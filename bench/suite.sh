#!/bin/bash
# Full-scale benchmark suite for the load-balancing work.
#
#   suite.sh [--wait]
#
# With --wait it blocks until no other lpjml process is running, so it can be
# queued behind a production run.  Timing on a shared machine is meaningless,
# so every step re-checks and refuses to start while something else is using
# the cores.
#
# Every run is derived from a real WHEP config.  The timing runs start from the
# 300-year spun-up restart, which is both faster and more representative than a
# fresh spinup: the cost of a cell depends on the vegetation standing on it.
set -u
export LC_ALL=C
cd "$(dirname "$0")/.." || exit 1
ROOT=$PWD
RUNS=/home/usuario/lpjml_perf_runs
REC=/home/usuario/lpjml_rerun_2026-08/recover_s5/config_recover.json
RECPREFIX=/home/usuario/lpjml_rerun_2026-08/recover_s5
COST=$RUNS/cost_measured.bin
LOG=$RUNS/suite.log

mkdir -p "$RUNS"
say() { echo "[$(date +%H:%M:%S)] $*" | tee -a "$LOG"; }

# Match the executable name, not the command line. Log paths, monitor loops and
# other sessions' shell wrappers all mention lpjml and would block this forever;
# one such monitor is what kept the first attempt waiting through an idle hour.
# Our own runs go through /tmp/lpjml_bench_bin, so they carry a different name.
lpjml_running() { pgrep -x lpjml >/dev/null || pgrep -x lpjml_timing >/dev/null; }

wait_for_idle() {
    while lpjml_running; do
        sleep 60
    done
    say "machine is idle"
    sleep 20   # let page cache and any stragglers settle
}

[ "${1:-}" = "--wait" ] && { say "waiting for the machine to go idle..."; wait_for_idle; }

# ---------------------------------------------------------------- P vs E core
# How much slower is an E-core than a P-core for this workload?  That ratio is
# what task_weights has to encode.  One rank pinned to CPU 0 (P-core) and one
# to CPU 16 (E-core), same work.
say "== P-core vs E-core ratio =="
"$ROOT/bench/mkconfig.py" --out "$RUNS/pe" --startgrid 30000 --endgrid 30499 \
    --river-routing false --nspinup 4 --firstyear 1901 --lastyear 1901 >/dev/null
pe_time() {
    cd "$RUNS/pe" && rm -f output/scenario_1/* restart/scenario_1/*
    local s e
    s=$(date +%s.%N)
    taskset -c "$1" "$ROOT/bin/lpjml" configurations/config_scenario_1.json >/dev/null 2>&1
    e=$(date +%s.%N)
    echo "$e - $s" | bc
    cd "$ROOT" || exit 1
}
T_P=$(pe_time 0)    # CPU 0  is a P-core thread
T_E=$(pe_time 16)   # CPU 16 is an E-core
PE_RATIO=$(echo "scale=3; $T_P / $T_E" | bc)
say "P-core ${T_P}s, E-core ${T_E}s -> an E-core is worth ${PE_RATIO} of a P-core"

# --------------------------------------------------------------- calibration
# One transient year from the spun-up restart, timing every cell, to produce
# the cost file the balanced runs use.
say "== calibration: measure per-cell cost (1 yr from restart, 24 ranks) =="
"$ROOT/bench/mkconfig.py" --out "$RUNS/cal" --base "$REC" --old-prefix "$RECPREFIX" \
    --lastyear 1901 --set "write_cellcost_filename=\"$COST\"" >/dev/null
"$ROOT/bench/run_bench.sh" "$RUNS/cal" "$ROOT/bin/lpjml" 24 | tee -a "$LOG"
"$ROOT/bench/cellcost.py" dump "$COST" | tee -a "$LOG"
"$ROOT/bench/cellcost.py" report "$COST" --ntask 24 | tee -a "$LOG"
"$ROOT/bench/cellcost.py" report "$COST" --ntask 30 | tee -a "$LOG"

# ------------------------------------------------ correctness with routing on
# The one gate still outstanding: does splitting by cost change anything when
# the tasks exchange river discharge every day?
say "== correctness: equal split vs cost split, river routing on =="
for n in rr_eq rr_bl; do
    "$ROOT/bench/mkconfig.py" --out "$RUNS/$n" --base "$REC" --old-prefix "$RECPREFIX" \
        --lastyear 1901 \
        --set 'write_restart=true' \
        --set "write_restart_filename=\"$RUNS/$n/restart/scenario_1/restart.lpj\"" >/dev/null
done
python3 - "$RUNS/rr_bl/configurations/config_scenario_1.json" "$COST" <<'PY'
import json, sys
p, cost = sys.argv[1], sys.argv[2]
c = json.load(open(p)); c["cellcost_filename"] = cost
json.dump(c, open(p, "w"), indent=1)
PY
"$ROOT/bench/run_bench.sh" "$RUNS/rr_eq" "$ROOT/bin/lpjml" 24 | tee -a "$LOG"
"$ROOT/bench/run_bench.sh" "$RUNS/rr_bl" "$ROOT/bin/lpjml" 24 | tee -a "$LOG"
say "-- comparison --"
"$ROOT/bench/cmp_runs.py" "$RUNS/rr_eq" "$RUNS/rr_bl" 2>&1 | tail -5 | tee -a "$LOG"
rm -f "$RUNS"/rr_eq/restart/scenario_1/* "$RUNS"/rr_bl/restart/scenario_1/*

# --------------------------------------------------------------- timing A / B
# Five transient years from restart, which is long enough for the fixed startup
# cost not to dominate.
say "== timing: equal split vs cost split =="
mkcfg() {  # name, extra json settings applied by python below
    "$ROOT/bench/mkconfig.py" --out "$RUNS/$1" --base "$REC" --old-prefix "$RECPREFIX" \
        --lastyear 1905 >/dev/null
}
addcost() {
    python3 - "$RUNS/$1/configurations/config_scenario_1.json" "$COST" "${2:-}" <<'PY'
import json, sys
p, cost = sys.argv[1], sys.argv[2]
w = sys.argv[3] if len(sys.argv) > 3 else ""
c = json.load(open(p)); c["cellcost_filename"] = cost
if w:
    c["task_weights"] = [float(x) for x in w.split(",")]
json.dump(c, open(p, "w"), indent=1)
PY
}

for n in t_eq24 t_bl24 t_eq30 t_bl30 t_bd24; do mkcfg "$n"; done
addcost t_bl24
addcost t_bl30
addcost t_bd24 "$(python3 -c "print(','.join(['1.0']*8+['${PE_RATIO}']*16))")"

"$ROOT/bench/run_bench.sh" "$RUNS/t_eq24" "$ROOT/bin/lpjml" 24 | tee -a "$LOG"
"$ROOT/bench/run_bench.sh" "$RUNS/t_bl24" "$ROOT/bin/lpjml" 24 | tee -a "$LOG"
"$ROOT/bench/run_bench.sh" "$RUNS/t_eq30" "$ROOT/bin/lpjml" 30 --use-hwthread-cpus | tee -a "$LOG"
"$ROOT/bench/run_bench.sh" "$RUNS/t_bl30" "$ROOT/bin/lpjml" 30 --use-hwthread-cpus | tee -a "$LOG"
"$ROOT/bench/run_bench.sh" "$RUNS/t_bd24" "$ROOT/bin/lpjml" 24 --bind-to core --map-by core | tee -a "$LOG"

say "-- outputs must still match --"
"$ROOT/bench/cmp_runs.py" "$RUNS/t_eq24" "$RUNS/t_bl24" 2>&1 | tail -3 | tee -a "$LOG"
"$ROOT/bench/cmp_runs.py" "$RUNS/t_eq24" "$RUNS/t_bd24" 2>&1 | tail -3 | tee -a "$LOG"

# ------------------------------------------------------------------- profile
say "== profile: where the time goes, and how far the tasks drift apart =="
mkcfg t_prof
addcost t_prof
"$ROOT/bench/run_bench.sh" "$RUNS/t_prof" "$ROOT/bin/lpjml_timing" 24 | tee -a "$LOG"
sed -n '/LPJmL performance summary/,$p' "$RUNS/t_prof/output/run.out" | tee -a "$LOG"

say "== suite complete =="
grep -h wall= "$LOG" | tail -20
