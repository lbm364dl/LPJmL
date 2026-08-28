# LPJmL performance work — measurement and correctness harness

The goal of this branch is to make the WHEP global run faster **without
changing a single number it produces**.  That is only tractable because of one
measured property of this model:

> LPJmL 6.1.1 on this fork is bit-reproducible, and its results do not depend
> on how the grid is split across MPI ranks.

Everything here exists to keep that property as an enforced gate rather than a
hope.

## The correctness gate

`cmp_runs.py A B` compares two run directories:

* every NetCDF output variable, value by value (not "close" — equal);
* every raw/text output, byte for byte;
* the restart file, which is the complete model state at the end of the run.

The restart header carries a wall-clock timestamp, a version string and a git
hash, so a run always differs there.  `cmp_runs.py` localises the differing
bytes and reports `IDENTICAL*` when they all fall inside the header, and
`DIFFERS` as soon as one byte of model state moves.  Runs are executed through
a fixed-length binary path (`/tmp/lpjml_bench_bin`) so that the header has the
same length for every binary being compared.

A change is accepted only if `cmp_runs.py` says `IDENTICAL`.  Any change that
cannot meet that bar is a scientific decision, not an optimisation, and is
handled separately.

### Established baselines (2026-08-28)

| comparison | result |
| --- | --- |
| same binary, same ranks, twice | 70/70 outputs identical, restart identical | 
| 8 ranks vs 24 ranks | 70/70 outputs identical, restart identical |

The second row is what makes load balancing a free change: the per-cell inflow
summation in `drain()` iterates `connect[i].index[j]` in an order fixed by the
routing network, and `pnet_setup()` only remaps the *values* of those indices
when the decomposition changes, never their order.

## Benchmarks

`mkconfig.py` derives every benchmark from the config that produced the real
11h29m global run, so the code path being measured is the one that matters:

    ./bench/mkconfig.py --out RUNDIR --startgrid 30000 --endgrid 35999 \
        --river-routing false --nspinup 10 --firstyear 1901 --lastyear 1902

`run_bench.sh RUNDIR BINARY RANKS [mpirun args...]` runs it and records wall
time, peak RSS and the mpirun arguments in `RUNDIR/bench.json`.

Cell subsets cannot use river routing — `initdrain()` fails when a cell drains
to a cell outside the subset — so routing benchmarks have to use the full grid.

## Measuring on this machine

Timing runs are only valid on an idle machine.  Check first:

    ps -eo pid,etime,pcpu,args | grep lpjml

An 11-hour production run and a benchmark on the same 24-core box produce
numbers that mean nothing.  Correctness comparisons are unaffected by load and
can be run at any time.

`perf` needs `kernel.perf_event_paranoid <= 1`; it is 4 here, so sampling
profiles require `sudo sysctl -w kernel.perf_event_paranoid=1` first.  Without
it, the built-in `-with_timing` build (`bin/lpjml_timing`) gives a per-function
breakdown with min/max/avg across ranks, which is also what exposes load
imbalance.
