#!/bin/bash
# Hunt the source of the run-to-run non-reproducibility of restart-based runs.
#
#   bench/find_nondeterminism.sh [reproduce|valgrind|asan]
#
# reproduce  run the 500-cell restart case N times and count distinct answers
# valgrind   run it once under valgrind and report only the LPJmL frames
# asan       run it once under AddressSanitizer
#
# The case is 500 cells, river routing off, one transient year from the 300-year
# spun-up restart: about five seconds, which is what makes this tractable.
# Everything larger reproduces the same behaviour and takes minutes.
#
# What is already known (see README.md): the state at the start of the year is
# identical between runs, the divergence appears during the first day's cell
# loop, and only in cells carrying five or more crop stands.  It is not
# uninitialised stack, and it does not need MPI or river routing.
set -u
export LC_ALL=C
cd "$(dirname "$0")/.." || exit 1
ROOT=$PWD
RUN=/home/usuario/lpjml_perf_runs/sub_a
REC=/home/usuario/lpjml_rerun_2026-08/recover_s5/config_recover.json
PRE=/home/usuario/lpjml_rerun_2026-08/recover_s5
MODE=${1:-reproduce}

[ -d "$RUN" ] || "$ROOT/bench/mkconfig.py" --out "$RUN" --base "$REC" --old-prefix "$PRE" \
    --startgrid 30000 --endgrid 30499 --river-routing false --lastyear 1901 >/dev/null
CFG=$RUN/configurations/config_scenario_1.json

build() {  # build a binary with the given extra flags, into $2
    local extra=$1 out=$2
    make clean >/dev/null 2>&1
    ./configure.sh -prefix "$ROOT" -inpath /home/usuario/WHEP/LPJmL_inputs \
        -noerror -DTRACE_DAILY >/dev/null 2>&1
    sed -i "s|^OPTFLAGS= -O3\$|OPTFLAGS= $extra|" Makefile.inc
    make -j"$(nproc)" >/dev/null 2>&1 || { echo "build failed" >&2; exit 1; }
    cp bin/lpjml "$out"
    make clean >/dev/null 2>&1
    ./configure.sh -prefix "$ROOT" -inpath /home/usuario/WHEP/LPJmL_inputs >/dev/null 2>&1
    make -j"$(nproc)" >/dev/null 2>&1
}

case $MODE in
reproduce)
    build "-O1 -g" /tmp/lpjml_nd
    cd "$RUN" || exit 1
    n=${2:-8}
    out=$(for i in $(seq "$n"); do
            rm -f output/scenario_1/*
            mpirun -np 4 /tmp/lpjml_nd "$CFG" 2>/dev/null |
                grep -a '^TRACE 1901   1 cells' | awk '{print $5}'
          done)
    echo "$out" | sort | uniq -c
    echo "-> $(echo "$out" | sort -u | wc -l) distinct results out of $n runs"
    ;;
valgrind)
    command -v valgrind >/dev/null || { echo "valgrind is not installed" >&2; exit 1; }
    build "-O1 -g" /tmp/lpjml_nd
    cd "$RUN" || exit 1
    rm -f output/scenario_1/*
    # one rank, so valgrind sees the model and not the MPI runtime
    valgrind --tool=memcheck --track-origins=yes --error-limit=no \
             --num-callers=25 --log-file=/tmp/vg.log \
             /tmp/lpjml_nd "$CFG" >/dev/null 2>&1
    echo "== uninitialised-value errors mentioning LPJmL source =="
    awk '/^==/{buf=buf"\n"$0} /Conditional jump|uninitialised|Invalid read|Invalid write/{flag=1}
         /^==[0-9]*== *$/{if(flag) print buf; buf=""; flag=0}' /tmp/vg.log |
        grep -B 2 -A 12 "LPJmL-perf/src" | head -80
    echo
    echo "full log: /tmp/vg.log"
    ;;
asan)
    build "-O1 -g -fsanitize=address" /tmp/lpjml_nd
    cd "$RUN" || exit 1
    rm -f output/scenario_1/*
    ASAN_OPTIONS=detect_leaks=0:halt_on_error=0 \
        mpirun -np 1 /tmp/lpjml_nd "$CFG" >/dev/null 2>/tmp/asan.log
    echo "AddressSanitizer errors: $(grep -ac 'ERROR: AddressSanitizer' /tmp/asan.log)"
    grep -a -A 12 'ERROR: AddressSanitizer' /tmp/asan.log | head -40
    ;;
*)
    echo "usage: $0 [reproduce|valgrind|asan]" >&2
    exit 1
    ;;
esac
