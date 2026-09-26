# DARTS-128 Key-Recovery PoC

## Principle

Step 9 of Algorithm 7 (Sign) of the [DARTS-128 specification](https://cdn.jsdelivr.net/gh/ngcc-dev/ngcc-harness@main/sign-08/sign-08-spec.pdf) compares $\mathrm{Comp}_d(\omega)$ with $\mathrm{Comp}_d(\omega + (1-b)\,c\alpha)$ and restarts if the two are not equal, where $b \in \{0,1\}$ is a branch bit sampled uniformly per attempt. For $b = 1$ the two arguments coincide, so this check never triggers a restart. For $b = 0$ it is a genuine test whose outcome depends on $\omega$, and hence on the secret key through $\omega = A\lfloor y \rceil$, so the $b = 0$ branch survives with a key-dependent probability. The two branches therefore survive with slightly different probabilities, and the joint distribution of the published signature values $(c, \lfloor z_1 \rceil)$ becomes tilted by a small multiple of the secret key $s$.

Averaging the statistic $C(c)^{\mathsf{T}} \lfloor z_1 \rceil$ over $N$ signatures yields $a s + \mu + \varepsilon$, where $a$ is a small positive amplitude, $\mu$ is a fixed key-dependent bias that is approximately *smooth* in the coefficient index and does not decrease with $N$, while $s$ is a sparse impulse train (a few hundred nonzero $\pm 1$ entries) and the noise decreases as $1/\sqrt{N}$. A cross-validated smoother combined with a ternary soft-EM separates the two components. A depth-bounded enumeration over the few most uncertain coefficients then corrects the remaining misclassifications: each candidate is checked against the exact public-key relation $e = \frac{q+1}{2}\,\beta - A s \pmod{q}$ with $\beta = (1, 0, \ldots, 0)$, and accepted only if the implied $e$ is ternary. Experiments with the reference implementation show that 20,000,000 to 30,000,000 signatures suffice to recover the full secret key of DARTS-128.

## Contents

| Path | Role |
|---|---|
| `build.sh` | Fetches the reference implementation into `./DARTS` if missing, copies the sources and `test_data/` into the DARTS128 directory, compiles `genkey`, `collect`, `dfs` |
| `genkey.c` | Key generation -> `pk.bin`, `sk.bin` |
| `collect.c` | Multithreaded signature collection -> `acc_a.bin`, `acc_b.bin` |
| `dfs.c` | Public-key oracle + depth-bounded enumeration engine (used by `recover.py`) |
| `recover.py` | Leakage detection, key-adaptive bias separation, sufficiency diagnosis, oracle correction, hash comparison |
| `test_data/` | Ready-made data set: default key pair + accumulators of 20,000,000 real signatures (skip collection for a quick test) |

## Steps

### 1. Prerequisites (Ubuntu/Debian)

```bash
sudo apt update && sudo apt install -y gcc python3-numpy unzip curl
```

### 2. Build

Unpack this package into an empty directory and run:

```bash
bash build.sh
```

If `./DARTS` is absent, `build.sh` downloads `DARTS.zip` from the [official site](https://www.niccs.org.cn/niccs/Proposal/Public-Key%20Cryptographic%20Algorithms/Round%201%20candidates/DARTS.zip), unpacks it, copies the PoC sources and `test_data/` into `DARTS/Implementations/Reference_Implementation/DARTS128/` and compiles the three C files there.

### 3. Generate the key pair

```bash
cd DARTS/Implementations/Reference_Implementation/DARTS128
./genkey
```

Expected output:

```
using reference-implementation default key (zero-init DRNG)
key generated: |s|_0 = 312 nonzeros, |e|_0 = 148 nonzeros
wrote pk.bin, sk.bin
```

The default key is deterministic, so `md5sum sk.bin` prints `559499f6a5ff3e28b7aea6da96f5371a` on every machine. 

To attack a different key: `POC_SEED=<any string> ./genkey`.

### 4. Collect signatures

```bash
./collect 20000000
```

- Uses all CPU cores automatically (override with `POC_THREADS=<n>`). Reference throughput on an 8-core cloud VM: ~950 sig/s, so 20M signatures take roughly 6 hours.
- Progress is checkpointed to `ckpt.bin` every 60 s. If the process is interrupted, rerun the same command to resume from the checkpoint.
- Collection is finished when it prints `wrote acc_a.bin, acc_b.bin`.

### 5. Recover the key

```bash
python3 recover.py --compare-sk
```

The script estimates the signature count from the noise between the odd/even accumulator halves, separates the key-dependent bias, searches for the exact key against the public-key oracle, and saves it as `recovered_sk.bin` (sk.bin format; the random seed K is not recoverable and is zeroed). With `--compare-sk` it additionally compares the result against `sk.bin` and prints the md5 of both.

Expected runtime: under a minute. Expected output:
```
signatures N ~ 19,983,582  (estimated from the accumulators)
leakage detection: z = +656.87  (|z| > 3 means significant)
stage 1 bias separation ... done
oracle search ... key found after 40,342 candidates (depth<=2)
recovered key written to recovered_sk.bin
comparing with sk.bin (md5 of bytes [0:1504]; last 32 bytes = signing seed K, excluded):
  sk.bin           = a08a9cc6dfc377a6f0fa1330ca63434a
  recovered_sk.bin = a08a9cc6dfc377a6f0fa1330ca63434a  <-- identical

POC PASSED: the DARTS-128 secret key was fully recovered from ~19,983,582 public signatures alone.
```

When the signature count is clearly too low, the script skips the futile search and instead prints the estimated requirement (rerun `./collect` with the larger count, which resumes from its checkpoint):

```
signature count insufficient: estimated need ~24,000,000 (collect ~19,000,000 more and re-run; ./collect resumes from its checkpoint)
```

The hash comparison excludes the last 32 bytes of `sk.bin` on purpose: they hold the signing-randomness seed `K`, which only feeds `H(K, H(pk, msg))`, is never checked during verification, and is irrelevant to forgery: with `(s0, s1, e)` any `K` signs valid messages.

## Quick test with the bundled data set

The `test_data/` directory contains `pk.bin`, `sk.bin` (the default key) and `acc_a.bin`, `acc_b.bin` (accumulators of 20,000,000 signatures collected with `collect` from the default key), so Step 4 can be skipped. `build.sh` copies it next to the compiled binaries, so right after building:

```bash
cd DARTS/Implementations/Reference_Implementation/DARTS128
python3 recover.py test_data --compare-sk
```

This finishes in about a minute with the `POC PASSED` line above.
