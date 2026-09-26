#!/usr/bin/env bash
# build.sh - fetches the reference implementation into ./DARTS if missing,
# copies the PoC sources into the DARTS128 directory and compiles
# genkey, collect and dfs.
set -euo pipefail

URL="https://www.niccs.org.cn/niccs/Proposal/Public-Key%20Cryptographic%20Algorithms/Round%201%20candidates/DARTS.zip"
SRC="DARTS/Implementations/Reference_Implementation/DARTS128"

if [ ! -d DARTS ]; then
    echo "[build] ./DARTS not found, preparing the reference implementation ..."
    if [ ! -f DARTS.zip ]; then
        echo "[build] downloading $URL"
        if command -v curl >/dev/null 2>&1; then
            curl -fL --retry 3 -o DARTS.zip "$URL"
        else
            wget -O DARTS.zip "$URL"
        fi
    else
        echo "[build] using existing ./DARTS.zip"
    fi
    mkdir -p DARTS
    unzip -q -o DARTS.zip -d DARTS
fi

if [ ! -d "$SRC" ]; then
    echo "[build] ERROR: $SRC not found after unpacking." >&2
    exit 1
fi

echo "[build] copying PoC sources into $SRC"
cp genkey.c collect.c dfs.c recover.py "$SRC"/
if [ -d test_data ]; then
    cp -r test_data "$SRC"/
fi
cd "$SRC"

IMPL_SRCS="packing.c polyvec.c poly.c ntt.c reduce.c sampler.c encoding.c \
polyfix.c polymat.c fft.c fixpoint.c symmetric.c drng.c auxfunc.c"

echo "[build] compiling genkey"
gcc -O3 -march=native -std=gnu99 -o genkey genkey.c sign.c $IMPL_SRCS -lm

echo "[build] compiling collect"
gcc -O3 -march=native -std=gnu99 -o collect collect.c sign.c $IMPL_SRCS -lm -lpthread

echo "[build] compiling dfs"
gcc -O3 -march=native -std=gnu99 -o dfs dfs.c $IMPL_SRCS -lm

echo
echo "[build] done."
echo
echo "Next steps:"
echo "    cd $SRC"
echo "    ./genkey"
echo "    ./collect 20000000"
echo "    python3 recover.py --compare-sk"
echo
echo "quick test with the bundled data set (no collection needed):"
echo "    cd $SRC"
echo "    python3 recover.py test_data --compare-sk"
