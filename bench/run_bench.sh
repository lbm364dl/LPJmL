#!/bin/bash
# Run one LPJmL benchmark and record wall time + a correctness fingerprint.
#
#   run_bench.sh <run_dir> <binary> <ranks> [mpirun extra args...]
#
# The run directory must already contain configurations/config_scenario_1.json
# (see mkconfig.py).  Everything is written inside the run directory so runs
# never interfere with each other.
set -u
export LC_ALL=C
RUN=$1; BIN=$2; NP=$3; shift 3
EXTRA=("$@")

CFG="$RUN/configurations/config_scenario_1.json"
[ -f "$CFG" ] || { echo "no config at $CFG" >&2; exit 1; }

mkdir -p "$RUN/output/scenario_1" "$RUN/restart/scenario_1"
# Deliberately does not delete anything.  LPJmL recreates its output files, and
# an unattended benchmark should never need a destructive action.  Use a fresh
# run directory per experiment (mkconfig.py --out) rather than reusing one.

# Exec through a fixed-length path so the argv string that LPJmL stamps into
# the restart-file header has the same length for every binary we compare.
STABLE=/tmp/lpjml_bench_bin
ln -sf "$(readlink -f "$BIN")" "$STABLE"

cd "$RUN" || exit 1
START=$(date +%s.%N)
nice -n 19 /usr/bin/time -v mpirun -np "$NP" "${EXTRA[@]}" "$STABLE" "$CFG" \
    >"$RUN/output/run.out" 2>"$RUN/output/run.err"
RC=$?
END=$(date +%s.%N)

WALL=$(echo "$END - $START" | bc)
MAXRSS=$(grep -a "Maximum resident set size" "$RUN/output/run.err" | tail -1 | tr -dc '0-9')
LPJWALL=$(grep -a "Total wall clock time" "$RUN/output/run.out" | tail -1 | sed 's/.*:\s*//')

cat >"$RUN/bench.json" <<JSON
{
 "run": "$(basename "$RUN")",
 "binary": "$BIN",
 "ranks": $NP,
 "mpirun_extra": "${EXTRA[*]}",
 "rc": $RC,
 "wall_s": $WALL,
 "lpjml_wall": "$LPJWALL",
 "max_rss_kb": ${MAXRSS:-0},
 "git": "$(git -C "$(dirname "$BIN")/.." rev-parse --short HEAD 2>/dev/null)"
}
JSON
printf '%-28s ranks=%-3s rc=%s wall=%.1fs  rss=%sMB  %s\n' \
    "$(basename "$RUN")" "$NP" "$RC" "$WALL" "$(( ${MAXRSS:-0} / 1024 ))" "${EXTRA[*]}"
exit $RC
