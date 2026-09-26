/* collect.c
 *
 * Description: Collects <num_sigs> signatures using all available CPU
 *   cores and accumulates the attack statistic into odd/even half-sums.
 *   Progress is checkpointed to ckpt.bin every 60 s; rerunning the same
 *   command after an interruption resumes from the checkpoint.
 *
 * Input:  sk.bin (from ./genkey); command-line argument <num_sigs>
 *
 * Output: acc_a.bin, acc_b.bin  (int64[L*N] half-sum accumulators)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include <unistd.h>
#include <signal.h>
#include <pthread.h>
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

/* caller must guarantee workers are all parked or all joined */
static int write_ckpt(worker_t *ws, long T, uint64_t nsig, uint64_t sksum) {
    FILE *f = fopen(CKPT_TMP, "wb");
    if (!f) return -1;
    fwrite(CKPT_MAGIC, 1, 8, f);
    fwrite(&nsig, sizeof nsig, 1, f);
    fwrite(&T, sizeof T, 1, f);
    fwrite(&sksum, sizeof sksum, 1, f);
    for (long i = 0; i < T; i++) {
        int64_t d = ws[i].done;
        fwrite(&d, sizeof d, 1, f);
        fwrite(ws[i].acc_a, sizeof(int64_t), L * N, f);
        fwrite(ws[i].acc_b, sizeof(int64_t), L * N, f);
    }
    fclose(f);
    return rename(CKPT_TMP, CKPT_FILE);       /* atomic replace */
}

/* returns 1 if a checkpoint was loaded, 0 for a fresh start */
static int load_ckpt(worker_t *ws, long T, uint64_t nsig, uint64_t sksum) {
    FILE *f = fopen(CKPT_FILE, "rb");
    if (!f) return 0;
    char magic[8]; uint64_t f_nsig, f_sksum; long f_T;
    if (fread(magic, 1, 8, f) != 8 || memcmp(magic, CKPT_MAGIC, 8) ||
        fread(&f_nsig, sizeof f_nsig, 1, f) != 1 || f_nsig != nsig ||
        fread(&f_T, sizeof f_T, 1, f) != 1 || f_T != T ||
        fread(&f_sksum, sizeof f_sksum, 1, f) != 1 || f_sksum != sksum) {
        fclose(f);
        fprintf(stderr, "warning: checkpoint does not match this run "
                        "(nsig/threads/key changed?), starting fresh\n");
        return 0;
    }
    uint64_t total = 0;
    for (long i = 0; i < T; i++) {
        int64_t d;
        if (fread(&d, sizeof d, 1, f) != 1 ||
            fread(ws[i].acc_a, sizeof(int64_t), L * N, f) != (size_t)(L * N) ||
            fread(ws[i].acc_b, sizeof(int64_t), L * N, f) != (size_t)(L * N)) {
            fclose(f);
            fprintf(stderr, "warning: checkpoint corrupted, starting fresh\n");
            /* reset progress only; begin/end/sk were set before the call */
            for (long k = 0; k < T; k++) {
                ws[k].done = 0;
                memset(ws[k].acc_a, 0, sizeof ws[k].acc_a);
                memset(ws[k].acc_b, 0, sizeof ws[k].acc_b);
            }
            return 0;
        }
        ws[i].done = d;
        total += (uint64_t)d;
    }
    fclose(f);
    fprintf(stderr, "resumed from checkpoint: %llu/%llu signatures done\n",
            (unsigned long long)total, (unsigned long long)nsig);
    return 1;
}

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "usage: %s <num_sigs>\n", argv[0]);
        return 1;
    }
    uint64_t nsig = strtoull(argv[1], NULL, 10);

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
    uint64_t chunk = nsig / (uint64_t)T, rem = nsig % (uint64_t)T, pos = 0;
    for (long i = 0; i < T; i++) {
        ws[i].begin = pos;
        pos += chunk + (i < (long)rem ? 1 : 0);
        ws[i].end = pos;
        ws[i].sk = sk;
    }

    int resumed = load_ckpt(ws, T, nsig, sk_checksum(sk));
    uint64_t done0 = 0;
    for (long i = 0; i < T; i++) done0 += (uint64_t)ws[i].done;

    signal(SIGINT, on_signal);
    signal(SIGTERM, on_signal);
    signal(SIGHUP, on_signal);

    struct timespec t0; clock_gettime(CLOCK_MONOTONIC, &t0);
    for (long i = 0; i < T; i++)
        if (pthread_create(&th[i], NULL, worker, &ws[i])) { fprintf(stderr, "pthread_create failed\n"); return 1; }

    /* progress monitor + periodic checkpointing */
    long finished = 0;
    time_t last_ckpt = time(NULL);
    while (finished < T && !stop_req) {
        usleep(5000000);
        uint64_t d = 0; finished = 0;
        for (long i = 0; i < T; i++) { d += (uint64_t)ws[i].done; if (ws[i].done == ws[i].end - ws[i].begin) finished++; }
        struct timespec t1; clock_gettime(CLOCK_MONOTONIC, &t1);
        double el = (t1.tv_sec - t0.tv_sec) + 1e-9 * (t1.tv_nsec - t0.tv_nsec);
        double rate = el > 0 ? (double)(d - done0) / el : 0;
        fprintf(stderr, "\r%llu/%llu sigs, %.0f sig/s, eta %.0f s   ",
                (unsigned long long)d, (unsigned long long)nsig, rate,
                rate > 0 ? (double)(nsig - d) / rate : 0.0);
        fflush(stderr);
        /* park all workers every 60 s and checkpoint atomically */
        if (finished < T && !stop_req && time(NULL) - last_ckpt >= CKPT_INTERVAL_SEC) {
            park_req = 1;
            while (parked < T && !stop_req) {
                long fin = 0;
                for (long i = 0; i < T; i++) if (ws[i].done == ws[i].end - ws[i].begin) fin++;
                if (fin == T) break;
                usleep(10000);
            }
            if (parked == T) {
                if (write_ckpt(ws, T, nsig, sk_checksum(sk)) == 0)
                    fprintf(stderr, "\n[checkpoint saved to %s]", CKPT_FILE);
            }
            park_req = 0;
            last_ckpt = time(NULL);
        }
    }
    for (long i = 0; i < T; i++) pthread_join(th[i], NULL);

    if (stop_req) {   /* interrupted: workers joined, state is consistent */
        uint64_t d = 0;
        for (long i = 0; i < T; i++) d += (uint64_t)ws[i].done;
        if (write_ckpt(ws, T, nsig, sk_checksum(sk)) == 0) {
            fprintf(stderr, "\ncheckpoint saved (%llu/%llu signatures). "
                            "Rerun the same command to resume.\n",
                    (unsigned long long)d, (unsigned long long)nsig);
        } else {
            fprintf(stderr, "\nfailed to save checkpoint!\n");
        }
        free(ws); free(th);
        return 3;
    }

    struct timespec t1; clock_gettime(CLOCK_MONOTONIC, &t1);
    double el = (t1.tv_sec - t0.tv_sec) + 1e-9 * (t1.tv_nsec - t0.tv_nsec);
    fprintf(stderr, "\ndone: %llu sigs in %.1fs (%.0f sig/s total, %.2f ms/sig/core)%s\n",
            (unsigned long long)nsig, el, nsig / el, 1000.0 * el / nsig * T,
            resumed ? " [resumed]" : "");

    /* merge per-thread accumulators and write out */
    static int64_t acc_a[L][N], acc_b[L][N];
    for (long i = 0; i < T; i++)
        for (int p = 0; p < L; p++)
            for (int j = 0; j < N; j++) {
                acc_a[p][j] += ws[i].acc_a[p][j];
                acc_b[p][j] += ws[i].acc_b[p][j];
            }
    f = fopen("acc_a.bin", "wb"); fwrite(acc_a, sizeof(int64_t), L * N, f); fclose(f);
    f = fopen("acc_b.bin", "wb"); fwrite(acc_b, sizeof(int64_t), L * N, f); fclose(f);
    unlink(CKPT_FILE);   /* completed: drop checkpoint */
    fprintf(stderr, "wrote acc_a.bin, acc_b.bin (nsig=%llu)\n",
            (unsigned long long)nsig);
    free(ws); free(th);
    return 0;
}
