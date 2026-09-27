/* dfs.c
 *
 * Description: Depth-bounded flip enumeration against the public-key
 *   oracle. A candidate s is valid iff r = HALF_Q*e0 - A*s (mod q) lies
 *   in {-1,0,1}^N (exact key-generation relation; A converted to the
 *   true coefficient domain). Candidates are evaluated incrementally
 *   in O(N) per node, with early termination; suspects are tried in
 *   file order (the caller sorts them by descending doubt).
 *
 * Input:  argv[1] = pk.bin; argv[2] = base candidate (int32[1024]);
 *         argv[3] = suspect positions (int32[m]); argv[4] = max depth;
 *         argv[5] = candidate budget
 *
 * Output: on hit, dfs_found.bin (int32[1024]) and dfs_found_e.bin
 *         (int32[512]); prints FOUND/NOTFOUND with the candidate count
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "params.h"
#include "packing.h"
#include "polyvec.h"
#include "polymat.h"
#include "reduce.h"

#define NN 512
#define LL 2
#define MONTINV 1535   /* 114369^{-1} mod 130817 */

static int32_t COL[LL][NN][NN];   /* negacyclic shifts of true-coefficient A (centered) */
static int32_t tvec[NN];
static int32_t base[LL * NN];
static int32_t susp[1024];
static int nsusp, maxdepth;
static uint64_t budget, tried;

static int32_t center_q(int64_t v) {
    v %= Q;
    if (v < 0) v += Q;
    if (v > Q / 2) v -= Q;
    return (int32_t)v;
}

static void build_columns(polyveck *A0, polyvecl A[K]) {
    polyveck A0c = *A0, A1c;
    for (int i = 0; i < K; i++) A1c.vec[i] = A[i].vec[1];
    polyveck_invntt_tomont(&A0c);
    polyveck_invntt_tomont(&A1c);
    int32_t Ac[LL][NN];
    for (int j = 0; j < NN; j++) {
        Ac[0][j] = center_q((int64_t)A0c.vec[0].coeffs[j] * MONTINV);
        Ac[1][j] = center_q((int64_t)A1c.vec[0].coeffs[j] * MONTINV);
    }
    for (int p = 0; p < LL; p++)
        for (int j = 0; j < NN; j++)
            for (int i = 0; i < NN; i++)
                COL[p][j][i] = (i >= j) ? Ac[p][i - j] : -Ac[p][i - j + NN];
}

static void compute_t(const int32_t *s) {
    memset(tvec, 0, sizeof tvec);
    for (int p = 0; p < LL; p++)
        for (int j = 0; j < NN; j++) {
            int32_t v = s[p * NN + j];
            if (!v) continue;
            for (int i = 0; i < NN; i++)
                tvec[i] = center_q((int64_t)tvec[i] + v * COL[p][j][i]);
        }
}

static int check(void) {
    for (int i = 0; i < NN; i++) {
        int32_t r = center_q((i == 0 ? HALF_Q : 0) - tvec[i]);
        if (r < -1 || r > 1) return 0;
    }
    return 1;
}

static int write_file(const char *path, const int32_t *buf, size_t n) {
    FILE *f = fopen(path, "wb");
    if (!f) return -1;
    int ok = fwrite(buf, sizeof(int32_t), n, f) == n;
    if (fclose(f) != 0) ok = 0;
    return ok ? 0 : -1;
}

/* returns 0 on success; on write failure the caller must not trust the files */
static int export_found(const int32_t *s) {
    int32_t earr[NN];
    for (int i = 0; i < NN; i++)
        earr[i] = center_q((i == 0 ? HALF_Q : 0) - tvec[i]);
    if (write_file("dfs_found.bin", s, LL * NN) ||
        write_file("dfs_found_e.bin", earr, NN)) {
        fprintf(stderr, "error: failed to write dfs_found*.bin\n");
        return -1;
    }
    return 0;
}

static int32_t cur[LL * NN];

static inline void apply(int p, int j, int d) {
    for (int i = 0; i < NN; i++)
        tvec[i] = center_q((int64_t)tvec[i] + d * COL[p][j][i]);
}

static int export_failed = 0;

static int dfs(int idx, int depth) {
    if (tried >= budget) return 0;
    if (depth == 0) {
        tried++;
        if (check()) {
            if (export_found(cur))
                export_failed = 1;
            return 1;
        }
        return 0;
    }
    for (int i = idx; i <= nsusp - depth; i++) {
        int j = susp[i];
        int p = j / NN, jj = j % NN;
        int32_t sj = cur[j];
        for (int v = -1; v <= 1; v++) {
            if (v == sj) continue;
            int d = v - sj;
            cur[j] = v;
            apply(p, jj, d);
            if (dfs(i + 1, depth - 1)) return 1;
            apply(p, jj, -d);
            cur[j] = sj;
        }
    }
    return 0;
}

int main(int argc, char **argv) {
    if (argc < 6) { fprintf(stderr, "usage: dfs <pk> <base> <suspects> <maxdepth> <budget>\n"); return 1; }
    uint8_t pk[CRYPTO_PUBLICKEYBYTES];
    FILE *f = fopen(argv[1], "rb");
    if (!f || fread(pk, 1, CRYPTO_PUBLICKEYBYTES, f) != CRYPTO_PUBLICKEYBYTES) return 1;
    fclose(f);
    polyveck A0; uint8_t seedA1[SEEDBYTES];
    unpack_pk(&A0, seedA1, pk);
    polyvecl A[K];
    polymatkl_expand(A, seedA1);
    for (int i = 0; i < K; i++) A[i].vec[0] = A0.vec[i];
    build_columns(&A0, A);

    f = fopen(argv[2], "rb");
    if (!f || fread(base, sizeof(int32_t), LL * NN, f) != LL * NN) return 1;
    fclose(f);
    f = fopen(argv[3], "rb");
    if (!f) return 1;
    nsusp = (int)(fread(susp, sizeof(int32_t), 1024, f));
    fclose(f);
    maxdepth = atoi(argv[4]);
    budget = strtoull(argv[5], NULL, 10);

    memcpy(cur, base, sizeof cur);
    compute_t(cur);
    tried = 1;
    if (check()) {
        if (export_found(cur))
            return 1;
        printf("FOUND depth=0 tried=1\n");
        return 0;
    }

    for (int d = 1; d <= maxdepth && tried < budget; d++) {
        if (dfs(0, d)) {
            if (export_failed) return 1;
            printf("FOUND depth<=%d tried=%llu\n", d, (unsigned long long)tried);
            return 0;
        }
        fprintf(stderr, "depth %d done, tried=%llu\n", d, (unsigned long long)tried);
    }
    printf("NOTFOUND tried=%llu\n", (unsigned long long)tried);
    return 2;
}
