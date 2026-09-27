/* genkey.c
 *
 * Description: Generates a DARTS-128 key pair with the reference
 *   implementation. Without POC_SEED the DRNG stays zero-initialized, so
 *   every run produces the same deterministic default key.
 *
 * Input:  none (optional environment variable POC_SEED selects another key)
 *
 * Output: pk.bin, sk.bin
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "params.h"
#include "sign.h"
#include "packing.h"
#include "polyvec.h"
#include "drng.h"

extern DRNG_ctx drng_algorithm;   /* defined in sign.c; zero-init => deterministic key */

static int write_file(const char *path, const uint8_t *buf, size_t len) {
    FILE *f = fopen(path, "wb");
    if (!f) return -1;
    int ok = fwrite(buf, 1, len, f) == len;
    if (fclose(f) != 0) ok = 0;
    return ok ? 0 : -1;
}

int main(void) {
    const char *env = getenv("POC_SEED");
    if (env && strlen(env) > 0) {
        uint8_t seed[32] = {0};
        size_t l = strlen(env);
        for (size_t i = 0; i < 32; i++) seed[i] = (uint8_t)env[i % l] ^ (uint8_t)(i * 0x9e);
        init_random_number(&drng_algorithm, seed, 32);
        fprintf(stderr, "POC_SEED set: generating a non-default key (deterministic per seed)\n");
    } else {
        fprintf(stderr, "using reference-implementation default key (zero-init DRNG)\n");
    }

    uint8_t pk[CRYPTO_PUBLICKEYBYTES], sk[CRYPTO_SECRETKEYBYTES];
    crypto_sign_keypair(pk, sk);

    if (write_file("pk.bin", pk, CRYPTO_PUBLICKEYBYTES) ||
        write_file("sk.bin", sk, CRYPTO_SECRETKEYBYTES)) {
        fprintf(stderr, "error: failed to write pk.bin/sk.bin\n");
        return 1;
    }

    /* report Hamming weights as a sanity check (truth stays inside sk.bin) */
    polyvecl A[K];
    poly s0;
    polyvecl_1 s1;
    polyveck e;
    uint8_t key[SEEDBYTES];
    unpack_sk(A, &s0, &s1, &e, key, sk);
    long w_s = 0, w_e = 0;
    for (int j = 0; j < N; j++) {
        w_s += (s0.coeffs[j] != 0);
        for (int p = 0; p < L - 1; p++) w_s += (s1.vec[p].coeffs[j] != 0);
        for (int p = 0; p < K; p++)     w_e += (e.vec[p].coeffs[j] != 0);
    }
    fprintf(stderr, "key generated: |s|_0 = %ld nonzeros, |e|_0 = %ld nonzeros\n", w_s, w_e);
    fprintf(stderr, "wrote pk.bin, sk.bin\n");
    return 0;
}
