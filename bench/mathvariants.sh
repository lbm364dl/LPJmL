#!/bin/bash
# Measure the compiler's floating-point options at one rank.
#
# The stock build targets generic x86-64, which means SSE2 and no FMA, on a
# chip that has AVX2 and FMA.  Three questions, in order of how much they are
# allowed to change the answer:
#
#   native   -march=native -ffp-contract=off : same arithmetic, wider registers
#   fma      -march=native                   : lets the compiler fuse a*b+c
#   fast     -march=native -ffast-math       : reassociation, no NaN/Inf care
#
# Each is timed and then gated against the stock build over the full restart
# state, so a change that is not identical still has to show how far it moved.
set -u
export LC_ALL=C
cd "$(dirname "$0")/.."
ROOT=$PWD
R=/home/usuario/lpjml_perf_runs
REC=/home/usuario/lpjml_rerun_2026-08/recover_s5/config_recover.json
PRE=/home/usuario/lpjml_rerun_2026-08/recover_s5
L=$R/mathvariants.log
: >"$L"
say(){ echo "[$(date +%H:%M:%S)] $*" | tee -a "$L"; }

exec 9>/tmp/lpjml_campaign.lock
flock -n 9 || { echo "another benchmark holds the lock" >&2; exit 1; }

mkgate() {  # name -> a run dir that writes a restart file
    local d=$R/gate_$1
    [ -d "$d" ] && return 0
    "$ROOT/bench/mkconfig.py" --out "$d" --base "$REC" --old-prefix "$PRE" \
        --startgrid 36000 --endgrid 36999 --river-routing false --lastyear 1903 >/dev/null
    python3 - "$d/configurations/config_scenario_1.json" <<'PY'
import json,sys
p=sys.argv[1]; c=json.load(open(p))
c["write_restart"]=True
c["write_restart_filename"]="restart/scenario_1/restart.lpj"
c["restart_year"]=1903
json.dump(c,open(p,"w"),indent=1)
PY
}

build() {  # extra flags -> /tmp/cs_$1
    make clean >/dev/null 2>&1
    ./configure.sh -prefix "$ROOT" -inpath /home/usuario/WHEP/LPJmL_inputs -noerror >/dev/null 2>&1
    [ -n "$2" ] && sed -i "s|^OPTFLAGS= -O3 -fno-math-errno\$|OPTFLAGS= -O3 -fno-math-errno $2|" Makefile.inc
    make -j"$(nproc)" >/tmp/build_$1.log 2>&1 || { say "  BUILD FAILED: $1"; return 1; }
    cp bin/lpjml /tmp/cs_$1
}

for v in "native:-march=native -ffp-contract=off" \
         "fma:-march=native" \
         "fast:-march=native -ffast-math"; do
    n=${v%%:*}; flags=${v#*:}
    say "== $n : $flags =="
    build "$n" "$flags" || continue
    mkgate "$n"
    for i in 1 2; do "$ROOT/bench/run_bench.sh" $R/gate_$n /tmp/cs_$n 1 | tee -a "$L"; done
    "$ROOT/bench/cmp_runs.py" $R/gate_lit $R/gate_$n 2>&1 | tail -6 | tee -a "$L"
done

# leave the tree as it ships
make clean >/dev/null 2>&1
./configure.sh -prefix "$ROOT" -inpath /home/usuario/WHEP/LPJmL_inputs -noerror >/dev/null 2>&1
make -j"$(nproc)" >/dev/null 2>&1
say "== math variants done =="
