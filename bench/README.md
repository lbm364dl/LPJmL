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

### Runs that start from a restart file did not reproduce themselves -- FIXED

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

**Cause**: `Cropdates` has seven members and the restart file carried six.
`fallow[2]`, the multicropping wait counter, was written by nobody and read by
nobody, so `freadcropdates()` left it as heap garbage and `sowingcft()` branched
on it on day one in every cell with cropland. Valgrind named it in one run.
Fixed by writing it to the restart and reading it back when present, so older
restart files still load.

After the fix, at full scale -- global grid, river routing, from the spun-up
restart -- these all hold:

| test | result |
| --- | --- |
| same configuration twice, 24 ranks | identical |
| equal split vs split by cost, 24 ranks | identical |
| 24 ranks vs 12 ranks | identical |

So the model is now reproducible, invariant to the number of tasks, and
invariant to how the grid is split. That is what makes byte-identical output a
usable gate for the performance work.

A canonical ordering of each cell's river inflows (`pnet_sortconnect`) was
written and then removed: rank-invariance holds without it, and it changed
results, so it was buying nothing.

### The load imbalance is real and large

From `bin/lpjml_timing`, one transient year, 24 ranks, equal cell counts:

    update_daily_cell   Tmin 5.57s   Tavg 46.37s   Tmax 89.16s

The heaviest task does sixteen times the work of the lightest. Splitting by
measured cell cost brings that to 26.61 / 45.72 / 75.82.

The residual is the machine: a P-core runs this workload in 21.8s where an
E-core takes 35.6s, so an E-core is worth 0.61 of a P-core. Ranks left on
E-cores set the pace no matter how the work is split, which is what
`task_weights` together with `--bind-to core --map-by core` is for.

## Measured result

Five transient years, full global grid, 24 ranks, from the spun-up restart.
Times are LPJmL's own "Wall clock time", so startup is excluded.

| configuration | simulation | speedup |
| --- | --- | --- |
| equal cell count, unpinned (current production setting) | 579 s | — |
| split by measured cell cost | 446 s | 1.30x |
| split by cost, pinned, weighted for P-cores vs E-cores | 380 s | 1.52x |

Then the source and compiler work (4.8% on its own) and the move to all 32
hardware threads with retuned weights:

| configuration | simulation | speedup |
| --- | --- | --- |
| equal cell count, 24 ranks unpinned | 546 s | — |
| split by cost, 32 threads pinned, weights retuned twice | **324 s** | **1.69x** |

0.0019 -> 0.0011 sec/cell/year.

With the full 71-output production set rather than the five-output recover
config, which is the load the real run carries:

| configuration | simulation | speedup |
| --- | --- | --- |
| equal cell count, 24 ranks unpinned | 579 s | — |
| split by cost, 32 threads pinned, weights retuned | **338 s** | **1.71x** |

70/70 output files identical. Against the 41,234 s production run that is
about **11.5 h -> 6.7 h**. Every step verified byte-identical.

The output set costs about 6% on a five-year run (579 s against 546 s for five
outputs), so the gather through the root task is real but not what limits it.

The cost file predicted 1.79x, and the gap is expected: 1.79x was the ceiling
for `update_daily_cell` alone, while river routing, output collection and the
monthly and annual steps do not shrink.

## Retuning the task weights

The weights say how fast each task is. An isolated core-speed measurement is
only a first guess: a P-core hyperthread pair, or sixteen E-cores sharing L2,
behave differently once every thread is busy. Measured here, the guess was
wrong by enough to matter -- P-threads actually spread over 0.86-1.00 and
E-cores over 0.57-0.79, a P/E ratio of 0.80 rather than the 0.92 an isolated
measurement suggests.

`cellcost.py retune` closes the loop. Run split by cost *and* writing a cost
file, then ask what each task's speed must have been given what it was handed
and how long it took:

    cellcost.py retune --split-by cost.bin --measured cost_insitu.bin \
        --ntask 32 --weights <the weights that run used>

It prints a new `task_weights` line. Two iterations converged here:

| iteration | slowest task / mean | wall |
| --- | --- | --- |
| first guess from isolated core speeds | 1.339 | 341 s |
| after one retune | 1.265 | 331 s |
| after two | **1.030** | **324 s** |

At 1.03 there is nothing left in load balancing; what remains is the
serialised part of the run.

## How to use it

Nothing in the model changes. Three optional config keys and one mpirun flag.

**Once**, to produce the cost file. Run one transient year from a spun-up
restart with the tasks on cores of the same kind, so that cell cost is not
confounded with core speed:

    "write_cellcost_filename": "/path/to/cost.bin"

    mpirun -np 8 --bind-to core --map-by core ./bin/lpjml config.json

Eight ranks land on the eight P-cores. It takes about three minutes.

`-fno-math-errno` is now in `config/Makefile.gcc` and `config/Makefile.mpich`,
so `configure.sh` picks it up.

**Then**, in every run:

    "cellcost_filename": "/path/to/cost.bin",
    "task_weights": [1,1,1,1,1,1,1,1, 0.61,0.61,0.61,0.61,0.61,0.61,0.61,0.61,
                     0.61,0.61,0.61,0.61,0.61,0.61,0.61,0.61]

    mpirun -np 32 --use-hwthread-cpus --bind-to hwthread --map-by hwthread \
        ./bin/lpjml config.json

Ranks 0-15 land two to a P-core and 16-31 one to an E-core. Start from the
isolated core speeds, then retune from a measured run (below) -- two iterations
is enough.

The weights are the relative speed of the core each rank is pinned to. With
`--bind-to core --map-by core` rank r goes to core r, and on this machine cores
0-7 are P-cores and 8-23 are E-cores; 0.61 is the measured ratio (a P-core runs
this workload in 21.8 s where an E-core takes 35.6 s). Re-measure on a
different machine. Both keys are optional and independent: with neither, the
split is exactly what it has always been.

Note that 24 pinned ranks beats the 30 unpinned ranks currently used. Going
past 24 needs `--use-hwthread-cpus` and a weight for every rank, and the second
thread of a P-core is worth much less than the first.

Recalibrate the cost file if the land use, the resolution or the machine
changes. A stale cost file costs speed and cannot change results.
