#!/bin/bash
# Build LPJmL with a profile, which is worth about 4.2% and is byte-identical.
#
#   bench/build_pgo.sh [training-config.json]
#
# Three steps: build instrumented, run once to collect a profile, rebuild
# against it with link-time optimisation as well.  The profile is not sensitive
# to which cells it was collected on -- one collected on dense agricultural
# cells is worth 4.2% there, 3.9% on a mixed block and 2.3% on a sparse one --
# so any representative year will do, and a short one is enough.  A wider
# training set, 40000 cells rather than 1000, is never worse and about half a
# percent better on cells unlike the training ones, so prefer a range spanning
# several latitudes if one is to hand.
#
# The profile lives in ./pgo and is reused until deleted.  Rebuild it after a
# change to the model, or the compiler will optimise for the old shape of the
# code; a stale profile costs performance and never accuracy.
set -eu
cd "$(dirname "$0")/.."
ROOT=$PWD
INPATH=${INPATH:-/home/usuario/WHEP/LPJmL_inputs}
TRAIN=${1:-}

if [ -z "$TRAIN" ]; then
    echo "usage: $0 <training-config.json>" >&2
    echo "  any representative configuration; one year is enough" >&2
    exit 1
fi
[ -f "$TRAIN" ] || { echo "no config at $TRAIN" >&2; exit 1; }

echo "== 1/3 instrumented build =="
make clean >/dev/null 2>&1
./configure.sh -prefix "$ROOT" -inpath "$INPATH" -noerror -pgo gen >/dev/null
make -j"$(nproc)" >/dev/null

echo "== 2/3 training run =="
( cd "$(dirname "$(dirname "$TRAIN")")" && "$ROOT/bin/lpjml" "$TRAIN" >/dev/null 2>&1 ) \
    || { echo "training run failed" >&2; exit 1; }
echo "   profile: $(ls pgo 2>/dev/null | wc -l) files in $ROOT/pgo"

echo "== 3/3 optimised build =="
make clean >/dev/null 2>&1
./configure.sh -prefix "$ROOT" -inpath "$INPATH" -noerror -lto -pgo use >/dev/null
make -j"$(nproc)" >/dev/null
echo "   done: $ROOT/bin/lpjml"
