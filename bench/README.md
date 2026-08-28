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

## The suite

`suite.sh --wait` queues the whole measurement behind whatever else is on the
machine and then runs, in order:

1. **P-core vs E-core.** One rank pinned to CPU 0, one to CPU 16, same work.
   The ratio is what `task_weights` has to encode.
2. **Calibration.** One transient year from the 300-year spun-up restart with
   `write_cellcost_filename` set, producing `cost_measured.bin`. Starting from
   a restart matters: what a cell costs depends on the vegetation standing on
   it, so cost measured from bare ground in spinup year 1 is the wrong profile.
3. **Correctness with river routing on.** Equal split vs cost split, both
   writing a restart, compared with `cmp_runs.py`. This is the one gate the
   small-scale tests could not cover, because a cell subset cannot use river
   routing -- `initdrain()` fails when a cell drains to a cell outside the
   subset.
4. **Timing A/B.** Five transient years, equal vs cost split, at 24 and 30
   ranks, plus a pinned-and-weighted run. Outputs are compared again after.
5. **Profile.** The same run under `bin/lpjml_timing`, whose per-function
   summary reports min/max/avg across tasks -- which is what shows how far the
   tasks still drift apart.

Predicted from the land-use proxy before any of this was measured: the current
equal-count split loads the heaviest of 30 tasks 1.50x above the mean, and
balancing recovers essentially all of it.

## Findings, 2026-08-28

### Runs that start from a restart file do not reproduce themselves

Two runs of the same configuration, same binary, same number of ranks, give
different results. Six runs of a 500-cell case produced five distinct answers.
This is **not** caused by anything on this branch: it reproduces on unmodified
`origin/master`.

What is known:

* A run that starts from bare ground (`restart: false`) is reproducible. Many
  runs, bit-identical. The trigger is reading a restart file.
* It is not river routing. It reproduces with `river_routing: false`.
* It is not MPI. It reproduces at one rank, just less often.
* The state at the start of the year is identical between runs -- carbon and
  nitrogen stocks, stand fractions, PFT counts, soil water, ice, and the soil
  heat grid all agree exactly. So the restart is read correctly.
* The divergence appears during the **first day's cell loop**, with identical
  precipitation, in cells carrying **five or more crop stands**. Cells with
  fewer never differ.
* It is not uninitialised stack: `-ftrivial-auto-var-init=pattern` does not fix
  it. `MALLOC_PERTURB_` results are inconclusive.
* Every field of `Pftcrop` and of `Irrigation` is restored by its reader.

Magnitude on a one-year run: about 1.6% of cells differ, some by tens of
percent locally; global means shift by 0.01-0.04%. Over 423 years the
divergence has 423 times as long to grow.

This matters beyond performance work: it means a production run cannot be
reproduced, and two runs that differ only in configuration cannot be compared
without knowing how much of the difference is noise.

`-DTRACE_DAILY` (`src/lpj/trace_day.c`) is the tool that localised it, and is
the tool for finishing the job. The next step is valgrind on the 500-cell
restart case, which is a five-second run.

### The load imbalance is real and large

From `bin/lpjml_timing`, one transient year, 24 ranks, equal cell counts:

    update_daily_cell   Tmin 5.57s   Tavg 46.37s   Tmax 89.16s

The heaviest task does sixteen times the work of the lightest. Splitting by
measured cell cost brings that to 26.61 / 45.72 / 75.82.

The residual is the machine: a P-core runs this workload in 21.8s where an
E-core takes 35.6s, so an E-core is worth 0.61 of a P-core. Ranks left on
E-cores set the pace no matter how the work is split, which is what
`task_weights` together with `--bind-to core --map-by core` is for.
