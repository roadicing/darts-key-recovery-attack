/* collect.c
 *
 * Description: Collects <num_sigs> signatures using all available CPU
 *   cores and accumulates the attack statistic into odd/even half-sums.
 *   Progress is checkpointed to ckpt.bin every 60 s. Without --load a
 *   run always starts over; with --load <file> the given checkpoint is
 *   used instead: the same count resumes it, a larger count first
 *   finishes the old collection and then appends fresh counters beyond
 *   it, so the final collection covers exactly [0, count) and the
 *   odd/even halves stay balanced; a smaller count is rejected. Only
 *   one collect process may run in a directory at a time (collect.lock).
 *
 * Input:  sk.bin (from ./genkey); arguments: <num_sigs> [--load <file>]
 *
 * Output: acc_a.bin, acc_b.bin  (int64[L*N] half-sum accumulators)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <errno.h>
#include <time.h>
#include <unistd.h>
#include <signal.h>
#include <pthread.h>
#include <fcntl.h>
#include <sys/file.h>
#include "params.h"
#include "sign.h"
#include "packing.h"
#include "polyvec.h"

#define CKPT_FILE "ckpt.bin"
#define CKPT_TMP  "ckpt.tmp"
#define CKPT_MAGIC "POCCKPT1"
#define CKPT_INTERVAL_SEC 60

typedef struct {
    uint64_t begin, end;            /* message counter range [begin, end) */
    const uint8_t *sk;
    volatile long done;             /* signatures completed (incl. resumed) */
    int64_t acc_a[L][N], acc_b[L][N];
} worker_t;

static volatile sig_atomic_t stop_req = 0;   /* signal received: stop now */
static volatile sig_atomic_t park_req = 0;   /* monitor asks workers to pause */
static volatile long parked = 0;

static void on_signal(int sig) {
    (void)sig;
    stop_req = 1;
    park_req = 0;                   /* release any parked workers */
}

static void *worker(void *arg) {
    worker_t *w = (worker_t *)arg;
    uint8_t sig[CRYPTO_SIGNATUREBYTES]; size_t siglen = 0;
    uint8_t msg[32];
    poly c; polyvecl low, high, z1; polyveck h;

    for (uint64_t ctr = w->begin + (uint64_t)w->done; ctr < w->end; ctr++) {
        if (stop_req) break;
        if (park_req && ((ctr & 0xff) == 0)) {  /* poll pause request every 256 sigs */
            __sync_fetch_and_add(&parked, 1);
            while (park_req && !stop_req) usleep(1000);
            __sync_fetch_and_sub(&parked, 1);
            if (stop_req) break;
        }
        memset(msg, 0, sizeof msg);
        memcpy(msg, &ctr, sizeof ctr);
        if (crypto_sign_signature(sig, &siglen, msg, sizeof msg, w->sk) != 0) { ctr--; continue; }
        if (unpack_sig(&c, &low, &high, &h, sig)) { fprintf(stderr, "unpack failed\n"); exit(1); }
        polyvecl_compose(&z1, &low, &high);   /* ~ round(z1), range (-64, Q) */

        int64_t (*acc)[N] = (ctr & 1) ? w->acc_b : w->acc_a;
        for (int p = 0; p < L; p++) {
            int32_t zc[N];
            for (int m = 0; m < N; m++) {
                int32_t v = z1.vec[p].coeffs[m];
                if (v > Q / 2) v -= Q;        /* center */
                zc[m] = v;
            }
            for (int l = 0; l < N; l++) {     /* negacyclic accumulation over supp(c) */
                if (!c.coeffs[l]) continue;
                for (int j = 0; j < N; j++) {
                    int m = j + l;
                    if (m >= N) acc[p][j] -= zc[m - N];
                    else        acc[p][j] += zc[m];
                }
            }
        }
        w->done++;
    }
    return NULL;
}

static uint64_t sk_checksum(const uint8_t *sk) {
    uint64_t s = 0;
    for (size_t i = 0; i < CRYPTO_SECRETKEYBYTES; i++) s = s * 131 + sk[i];
    return s;
}

/* split the counter range [lo, hi) over the T workers */
static void partition(worker_t *ws, long T, uint64_t lo, uint64_t hi) {
    uint64_t total = hi - lo, pos = lo;
    uint64_t chunk = total / (uint64_t)T, rem = total % (uint64_t)T;
    for (long i = 0; i < T; i++) {
        ws[i].begin = pos;
        pos += chunk + (i < (long)rem ? 1 : 0);
        ws[i].end = pos;
        ws[i].done = 0;
    }
}

static int write_vec_file(const char *path, const int64_t *v) {
    FILE *f = fopen(path, "wb");
    if (!f) return -1;
    int ok = fwrite(v, sizeof(int64_t), (size_t)(L * N), f) == (size_t)(L * N);
    if (fclose(f) != 0) ok = 0;
    return ok ? 0 : -1;
}

/* merge per-worker accumulators and write acc_a.bin/acc_b.bin */
static int export_acc(const worker_t *ws, long T) {
    static int64_t acc_a[L][N], acc_b[L][N];
    for (long i = 0; i < T; i++)
        for (int p = 0; p < L; p++)
            for (int j = 0; j < N; j++) {
                acc_a[p][j] += ws[i].acc_a[p][j];
                acc_b[p][j] += ws[i].acc_b[p][j];
            }
    if (write_vec_file("acc_a.bin", (const int64_t *)acc_a) ||
        write_vec_file("acc_b.bin", (const int64_t *)acc_b)) {
        fprintf(stderr, "error: failed to write acc_a.bin/acc_b.bin\n");
        return -1;
    }
    return 0;
}

/* caller must guarantee every worker is parked, finished, or joined */
static int write_ckpt(worker_t *ws, long T, uint64_t nsig, uint64_t sksum) {
    FILE *f = fopen(CKPT_TMP, "wb");
    if (!f) return -1;
    int ok = fwrite(CKPT_MAGIC, 1, 8, f) == 8 &&
             fwrite(&nsig, sizeof nsig, 1, f) == 1 &&
             fwrite(&T, sizeof T, 1, f) == 1 &&
             fwrite(&sksum, sizeof sksum, 1, f) == 1;
    for (long i = 0; i < T && ok; i++) {
        int64_t d = ws[i].done;
        ok = fwrite(&ws[i].begin, sizeof(uint64_t), 1, f) == 1 &&
             fwrite(&ws[i].end, sizeof(uint64_t), 1, f) == 1 &&
             fwrite(&d, sizeof d, 1, f) == 1 &&
             fwrite(ws[i].acc_a, sizeof(int64_t), L * N, f) == (size_t)(L * N) &&
             fwrite(ws[i].acc_b, sizeof(int64_t), L * N, f) == (size_t)(L * N);
    }
    if (fclose(f) != 0) ok = 0;
    if (!ok) {           /* keep the previous checkpoint, drop the partial */
        unlink(CKPT_TMP);
        return -1;
    }
    return rename(CKPT_TMP, CKPT_FILE);       /* atomic replace */
}

/* loads the checkpoint at <path> for a --load run. Returns 1 for a
   same-count resume and 2 when extending to a larger count; -1 means
   the checkpoint cannot be used. For an extension the checkpoint's
   ranges are kept as-is: the run first finishes the old collection,
   then appends fresh counters beyond it (two-phase; *old_nsig receives
   the checkpoint's bound count), so the final collection covers exactly
   the contiguous counter range [0, nsig). */
static int load_ckpt(worker_t *ws, long T, uint64_t nsig, uint64_t sksum,
                     const char *path, uint64_t *old_nsig) {
    FILE *f = fopen(path, "rb");
    if (!f) {
        fprintf(stderr, "error: cannot open checkpoint %s\n", path);
        return -1;
    }
    char magic[8]; uint64_t f_nsig, f_sksum; long f_T;
    int hdr_ok = fread(magic, 1, 8, f) == 8 && !memcmp(magic, CKPT_MAGIC, 8) &&
                 fread(&f_nsig, sizeof f_nsig, 1, f) == 1 &&
                 fread(&f_T, sizeof f_T, 1, f) == 1 &&
                 fread(&f_sksum, sizeof f_sksum, 1, f) == 1;
    if (!hdr_ok) {
        fclose(f);
        fprintf(stderr, "error: %s is not a valid checkpoint file\n", path);
        return -1;
    }
    if (f_T != T || f_sksum != sksum) {
        fclose(f);
        fprintf(stderr, "error: checkpoint belongs to an incompatible run "
                        "(different key or thread count)\n");
        return -1;
    }
    if (f_nsig > nsig) {
        fclose(f);
        fprintf(stderr, "error: checkpoint is bound to %llu signatures; "
                        "the count must be at least that (got %llu)\n",
                (unsigned long long)f_nsig, (unsigned long long)nsig);
        return -1;
    }
    int complete = 1;
    uint64_t total = 0, range_total = 0;
    for (long i = 0; i < T; i++) {
        uint64_t b, e; int64_t d;
        if (fread(&b, sizeof b, 1, f) != 1 || fread(&e, sizeof e, 1, f) != 1 ||
            fread(&d, sizeof d, 1, f) != 1 ||
            fread(ws[i].acc_a, sizeof(int64_t), L * N, f) != (size_t)(L * N) ||
            fread(ws[i].acc_b, sizeof(int64_t), L * N, f) != (size_t)(L * N)) {
            fclose(f);
            fprintf(stderr, "error: checkpoint %s is truncated or corrupted\n",
                    path);
            return -1;
        }
        ws[i].begin = b;
        ws[i].end = e;
        ws[i].done = d;
        total += (uint64_t)d;
        range_total += e - b;
        if (d != (int64_t)(e - b)) complete = 0;
    }
    fclose(f);
    /* the checkpoint's sums cover f_nsig - range_total counters from
       earlier runs plus the done prefixes of the current ranges */
    if (total > range_total || range_total > f_nsig) {
        fprintf(stderr, "error: checkpoint %s is truncated or corrupted\n",
                path);
        return -1;
    }
    *old_nsig = f_nsig;
    if (f_nsig == nsig) {
        /* counters below nsig - range_total were summed by earlier runs */
        fprintf(stderr, "resumed from checkpoint: %llu/%llu signatures done\n",
                (unsigned long long)(nsig - range_total + total),
                (unsigned long long)nsig);
        return 1;
    }
    /* larger count: keep the checkpoint's ranges and finish the old
       collection first; the monitor then assigns the fresh range
       [f_nsig, nsig). The final collection covers exactly [0, nsig), so
       the odd/even halves stay balanced (counters are never reused).
       Each phase announces the count its checkpoints are bound to. */
    if (complete)
        fprintf(stderr, "extending a completed collection: %llu -> %llu "
                        "signatures\n",
                (unsigned long long)f_nsig, (unsigned long long)nsig);
    else
        fprintf(stderr, "extending a partial collection: %llu -> %llu "
                        "signatures\n",
                (unsigned long long)(f_nsig - range_total + total),
                (unsigned long long)f_nsig);
    return 2;
}

int main(int argc, char **argv) {
    const char *num = NULL, *load_path = NULL;
    int bad = 0;
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--load")) {
            if (i + 1 >= argc || load_path) bad = 1;
            else load_path = argv[++i];
        } else if (!num) {
            num = argv[i];
        } else {
            bad = 1;
        }
    }
    if (!num || bad) {
        fprintf(stderr, "usage: %s <num_sigs> [--load <checkpoint>]\n", argv[0]);
        return 1;
    }
    char *endp = NULL;
    errno = 0;
    uint64_t nsig = strtoull(num, &endp, 10);
    /* strtoull accepts a leading '-' and wraps, so reject non-digits first */
    if (num[0] < '0' || num[0] > '9' || errno || !endp || *endp ||
        endp == num || nsig == 0) {
        fprintf(stderr, "invalid signature count: %s\n", num);
        return 1;
    }

    /* one collect per directory; the lock is held until process exit */
    int lockfd = open("collect.lock", O_RDWR | O_CREAT, 0644);
    if (lockfd < 0) {
        fprintf(stderr, "error: cannot open collect.lock: %s\n",
                strerror(errno));
        return 1;
    }
    if (flock(lockfd, LOCK_EX | LOCK_NB) != 0) {
        fprintf(stderr, "error: another collect process is already running "
                        "in this directory\n");
        return 1;
    }

    uint8_t sk[CRYPTO_SECRETKEYBYTES];
    FILE *f = fopen("sk.bin", "rb");
    if (!f || fread(sk, 1, CRYPTO_SECRETKEYBYTES, f) != CRYPTO_SECRETKEYBYTES) {
        fprintf(stderr, "sk.bin missing/bad - run ./genkey first\n"); return 1;
    }
    if (f) fclose(f);

    long T = sysconf(_SC_NPROCESSORS_ONLN);
    const char *env = getenv("POC_THREADS");
    if (env) T = atol(env);
    if (T < 1) T = 1;
    fprintf(stderr, "using %ld threads for %llu signatures\n", T, (unsigned long long)nsig);

    worker_t *ws = calloc(T, sizeof(worker_t));
    pthread_t *th = malloc(sizeof(pthread_t) * T);
    for (long i = 0; i < T; i++) ws[i].sk = sk;
    partition(ws, T, 0, nsig);

    /* state: 0 fresh, 1 resumed, 2 extension finishing the old collection */
    int state = 0;
    uint64_t old_nsig = 0;
    if (load_path) {
        state = load_ckpt(ws, T, nsig, sk_checksum(sk), load_path,
                          &old_nsig);
        if (state < 0) {
            free(ws); free(th);
            return 1;
        }
    } else if (access(CKPT_FILE, F_OK) == 0) {
        fprintf(stderr, "ignoring existing ckpt.bin (use --load ckpt.bin to "
                        "resume it)\n");
    }
    uint64_t done0 = 0, range_total = 0;
    for (long i = 0; i < T; i++) {
        done0 += (uint64_t)ws[i].done;
        range_total += ws[i].end - ws[i].begin;
    }
    /* an extension runs in two phases: first the old collection is
       finished (checkpoints stay bound to the old count), then the fresh
       range beyond it is collected and the checkpoint is bound to nsig */
    uint64_t bound = (state == 2) ? old_nsig : nsig;
    uint64_t phase1_new = 0;    /* signatures collected in phase 1 */
    /* progress display always shows the true summed total: counters
       below bound - range_total were summed by earlier runs */
    uint64_t disp_base = bound - range_total;
    uint64_t disp_target = nsig;
    /* extension of an already-complete collection: switch to the fresh
       range immediately (idempotent - the on-disk checkpoint still
       describes the old collection until the next save) */
    if (state == 2 && done0 == range_total) {
        partition(ws, T, old_nsig, nsig);
        done0 = 0;
        range_total = nsig - old_nsig;
        bound = nsig;
        disp_base = old_nsig;
        state = 1;
    }

    signal(SIGINT, on_signal);
    signal(SIGTERM, on_signal);
    signal(SIGHUP, on_signal);

    struct timespec t0; clock_gettime(CLOCK_MONOTONIC, &t0);
    int threads_up = 0;
    for (long i = 0; i < T; i++)
        if (pthread_create(&th[i], NULL, worker, &ws[i])) { fprintf(stderr, "pthread_create failed\n"); return 1; }
    threads_up = 1;

    /* progress monitor + periodic checkpointing */
    time_t last_ckpt = time(NULL);
    for (;;) {
        usleep(5000000);
        uint64_t d = 0;
        long finished = 0;
        for (long i = 0; i < T; i++) {
            d += (uint64_t)ws[i].done;
            if ((uint64_t)ws[i].done == ws[i].end - ws[i].begin) finished++;
        }
        struct timespec t1; clock_gettime(CLOCK_MONOTONIC, &t1);
        double el = (t1.tv_sec - t0.tv_sec) + 1e-9 * (t1.tv_nsec - t0.tv_nsec);
        double rate = el > 0 ? (double)(phase1_new + d - done0) / el : 0;
        fprintf(stderr, "\r%llu/%llu sigs, %.0f sig/s, eta %.0f s   ",
                (unsigned long long)(disp_base + d),
                (unsigned long long)disp_target, rate,
                rate > 0 ? (double)(disp_target - disp_base - d) / rate : 0.0);
        fflush(stderr);
        if (finished == T) {
            if (state != 2 || stop_req) break;  /* collection complete */
            /* phase 1 of an extension complete: join the workers before
               replacing their ranges - a worker may still be between its
               last done++ and the loop exit when done reaches the range
               end, and only a join guarantees nobody is in flight */
            for (long i = 0; i < T; i++) pthread_join(th[i], NULL);
            threads_up = 0;
            if (stop_req) break;
            phase1_new = range_total - done0;
            fprintf(stderr, "\nextending a completed collection: %llu -> "
                            "%llu signatures\n",
                    (unsigned long long)bound, (unsigned long long)nsig);
            partition(ws, T, bound, nsig);
            bound = nsig;
            done0 = 0;
            range_total = nsig - old_nsig;
            disp_base = old_nsig;
            state = 1;
            /* bind the checkpoint to the new count before starting the
               phase-2 workers, so an interrupt there resumes correctly */
            if (write_ckpt(ws, T, bound, sk_checksum(sk)) == 0)
                fprintf(stderr, "[checkpoint saved to %s]", CKPT_FILE);
            for (long i = 0; i < T; i++)
                if (pthread_create(&th[i], NULL, worker, &ws[i])) { fprintf(stderr, "pthread_create failed\n"); return 1; }
            threads_up = 1;
            last_ckpt = time(NULL);
            continue;
        }
        if (stop_req) break;
        /* park all workers every 60 s and checkpoint atomically */
        if (time(NULL) - last_ckpt >= CKPT_INTERVAL_SEC) {
            park_req = 1;
            /* wait until every worker is parked or finished; a finished
               worker never parks, so parked == T alone cannot be the goal */
            for (;;) {
                long fin = 0;
                for (long i = 0; i < T; i++)
                    if ((uint64_t)ws[i].done == ws[i].end - ws[i].begin) fin++;
                if (parked + fin == T || stop_req) break;
                usleep(10000);
            }
            if (!stop_req) {
                if (write_ckpt(ws, T, bound, sk_checksum(sk)) == 0)
                    fprintf(stderr, "\n[checkpoint saved to %s]", CKPT_FILE);
            }
            park_req = 0;
            last_ckpt = time(NULL);
        }
    }
    if (threads_up)
        for (long i = 0; i < T; i++) pthread_join(th[i], NULL);

    if (stop_req) {   /* interrupted: workers joined, state is consistent */
        uint64_t d = 0;
        for (long i = 0; i < T; i++) d += (uint64_t)ws[i].done;
        if (write_ckpt(ws, T, bound, sk_checksum(sk)) == 0) {
            fprintf(stderr, "\ncheckpoint saved (%llu/%llu signatures). "
                            "Resume with './collect %llu --load %s'.\n",
                    (unsigned long long)(disp_base + d),
                    (unsigned long long)disp_target,
                    (unsigned long long)nsig, CKPT_FILE);
        } else {
            fprintf(stderr, "\nfailed to save checkpoint!\n");
        }
        free(ws); free(th);
        return 3;
    }

    struct timespec t1; clock_gettime(CLOCK_MONOTONIC, &t1);
    double el = (t1.tv_sec - t0.tv_sec) + 1e-9 * (t1.tv_nsec - t0.tv_nsec);
    /* signatures collected by this run */
    uint64_t new_sigs = phase1_new + range_total - done0;
    if (new_sigs > 0 && el > 0)
        fprintf(stderr, "\ndone: %llu sigs in %.1fs (%.0f sig/s, %.2f ms/sig/core)%s\n",
                (unsigned long long)disp_target, el, new_sigs / el,
                1000.0 * el / new_sigs * T, state ? " [resumed]" : "");
    else
        fprintf(stderr, "\ndone: %llu sigs%s\n",
                (unsigned long long)disp_target, state ? " [resumed]" : "");

    /* keep the final checkpoint: a later --load run with a larger count
       extends it; saved before the export so that an export failure can
       still be recovered by simply rerunning with --load */
    if (write_ckpt(ws, T, bound, sk_checksum(sk)) != 0)
        fprintf(stderr, "warning: final checkpoint not saved; a later "
                        "larger-count run cannot extend this collection\n");

    if (export_acc(ws, T) != 0) {
        free(ws); free(th);
        return 1;
    }
    fprintf(stderr, "wrote acc_a.bin, acc_b.bin (nsig=%llu)\n",
            (unsigned long long)disp_target);
    free(ws); free(th);
    return 0;
}
