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

`perf` needs `kernel.perf_event_paranoid <= 1`; it has been set to 1 here, so
sampling profiles work directly.  Without
it, the built-in `-with_timing` build (`bin/lpjml_timing`) gives a per-function
breakdown with min/max/avg across ranks, which is also what exposes load
imbalance.

### The chip is power-limited, and that invalidates multi-rank A/B timing

The package runs under a sustained limit of 111 W averaged over an 80-second
window, with 253 W available in short bursts:

    /sys/class/powercap/intel-rapl:0/constraint_0_power_limit_uw -> 111 W / 80 s
    /sys/class/powercap/intel-rapl:0/constraint_1_power_limit_uw -> 253 W / 2 ms

An all-core run therefore starts fast on a full budget and settles to roughly
2.1 GHz once the moving average catches up, against 3.7 GHz when the budget has
had time to recover.  How fast a run goes depends on how long the machine idled
before it.  The same binary and configuration on 32 ranks produced 324, 328,
331, 332, 387 and 565 seconds -- a 1.7x spread with nothing changed.  It is not
thermal: the 565-second run was at 71 C and 3.7 GHz while an earlier 332-second
run was at 52 C and 2.1 GHz, hotter and faster at once.

This is why the earlier `nchunk` comparisons contradicted each other.  With a
fixed A/B ordering, whichever variant ran first got the recovered budget, and
the ordering bias was larger than the effect being measured.

The consequence for this harness: **do not A/B compare 32-rank runs.**  A single
rank draws about 15 W, far below the ceiling, and repeats to about 0.1 s.  Every
change to per-cell work -- which is all of the compute-side work below -- is
measured perfectly well, and far more honestly, on one rank.  Multi-rank runs
remain the right tool for load-balance questions, where the effect is large, but
each arm needs its ordering reversed and repeated.

## Findings, 2026-08-29: where the time actually goes

A sampling profile of a single rank, which had not been done before, moved the
work off guesswork.  The stock build spends about an eighth of its time inside
`pow()`:

    12.6%  littersom_nomethane
    12.4%  pow  (__ieee754_pow_fma 9.5% + pow@GLIBC 2.9%)
     3.4%  photosynthesis
     1.9%  exp
     ~8%   dynamic linker and MPI startup -- an artifact of a 62 s run

Three call sites accounted for most of the `pow()` time, and all three were
asking for the same answer repeatedly.  Each change below is exact: it removes a
recomputation, and never alters an arithmetic expression.

  * `photosynthesis()` computes `ko`, `kc` and `tau` as `pow()` of a fixed base
    against `(temp-25)*0.1`, and they feed only `fac` and `gammastar`.  Nothing
    there depends on anything but the temperature -- while `water_stressed()`
    holds the temperature fixed and calls the function up to fifty times as it
    bisects for lambda, and every stand in a cell sees the same daily
    temperature.  Remembering the last temperature returns the same bits.

  * The litter loop in `littersom_nomethane()` calls `pow(q10_wood, e)` twice
    per PFT with `e` fixed for the whole loop, and the PFT table holds only five
    distinct `q10_wood` values -- of which 1.0 covers every crop.  Answering
    each distinct base once per loop is exact.

  * `petpar()` recomputed the solar geometry, the longwave correction and the
    vapour-pressure slope for every stand, though albedo is its only per-stand
    input.  Split into `petpar_cell()` and `petpar_stand()`.

Measured at one rank on 1000 cells over three years, against the full 59 MB
restart state and all five output files:

    stock                                     62.8 s
    + petpar split                            62.0 s   -1.3%   identical
    + photosynthesis temperature cache        59.4 s   -5.4%   identical
    + littersom q10 table                     58.2 s   -7.3%   identical

What is left of `pow()` after those three is spread thin -- `f_wfps` ~5%,
`calc_soil_thermal_props` ~2.5%, `volatilization` ~2.2%, `infil_perc` ~1.9%,
`pedotransfer` ~1.5% -- with no exact reduction available.  `f_wfps` raises a
varying base to a per-soil-type real exponent, and `pedotransfer` looks like a
constant-per-cell function but depends on soil carbon, which evolves daily.

### Confirmed on the production configuration

The measurements above are 1000 cells with river routing off, because routing
needs the whole grid -- the irrigation neighbours of a subset point outside it,
and LPJmL refuses to build the network.  The changes were therefore checked
again on the configuration that is actually run: the full global grid, river
routing on, one year, 24 ranks, master against the optimised build, in both
orders so the power budget cannot favour either.

    master first     master 107.4 s   optimised  92.9 s    -13.5%
    optimised first  master 108.7 s   optimised 101.0 s     -7.1%
    ---------------------------------------------------------------
    mean             master 108.1 s   optimised  97.0 s    -10.3%

Byte-identical in both arms, all five output files.  The spread between the two
arms is the power-budget noise described above, about +/-3%, which is why both
orders are needed to state a number at all.  Roughly an hour off an 11.5-hour
run.

Note that 24 is the rank ceiling here: mpirun counts the 8 P-cores and 16
E-cores as 24 slots and refuses -np 32 without oversubscribing.

### Finding the redundant work instead of guessing at it

The three call sites in the first round were found by reading a profile and
reasoning about them.  That does not scale, and the two guesses that followed
were both wrong: a `pow(10,0)` guard that duplicated a fast path glibc already
has (21.7 million calls, all skipped, no measurable gain), and a `pow(1,y)`
guard whose condition never fires at all.

`bench/powprof.h` replaces the guessing.  It wraps every `pow()`/`exp()` call
site and records how often it is handed the argument it was handed last time.
A site with a high repeat rate is one a one-entry memo collapses exactly.

    ./configure.sh ... && sed -i 's|^OPTFLAGS= .*|& -DPOW_PROFILE -include '$PWD'/bench/powprof.h|' Makefile.inc
    make clean && make -j$(nproc)          # a clean build: it is a header change

Each translation unit keeps its own table and prints at exit, so nothing has to
be linked in.  **Build clean**, or only the files that happen to be newer than
their objects get instrumented -- which is how the first run of it came back
claiming every call in the model was in `infil_perc.c`.

It found five sites in the litter loop at 105 million calls each, repeating
between 57% and 98% of the time, and a `pow()` of a compile-time constant
sitting inside a loop in two places.  1579 million libm calls became 1176
million; `pow` fell from 9.2% of the run to 2.2% and `exp` left the top twenty.

A caution the tool cannot give you: a memo only pays when what it skips costs
more than the lookup.  An eight-entry cache would have raised the hit rate on
`temp_stress` from 37% to 80%, but a linear scan of eight doubles costs about
as much as the `exp()` it avoids, so it would have bought nothing.  The wins
here are the ones where the skipped work is a whole `pow()` or the lookup is
free because the value was hoisted out of a loop entirely.

### The SAFE checks cost about 4%

`-DSAFE` is on by default in every makefile template.  Its blocks are pure
checks -- they test for states that should not arise, and call `fail()` if one
does -- and there are 78 of them, 19 in `littersom_nomethane` and `infil_perc`
alone.

    with -DSAFE     56.0 s
    without         53.9 s     -3.7%     byte-identical

Byte-identical over the full restart state, which is what one would expect:
while no check fires, removing them cannot change an answer.  `./configure.sh
-nosafe` now turns them off.

It is off by default and should stay that way for anything exploratory.  The
checks are the difference between a run that stops and tells you the soil NO3
went negative, and a run that finishes and hands you the output anyway.  For a
configuration that has already been run through successfully, 4% for a
consistency check you have already passed is a reasonable trade; for a new one
it is not.

### The compiler's floating-point options

Measured the same way, at one rank, against the byte-identical build at 58.2 s.
The stock build targets generic x86-64 -- SSE2, no FMA -- on a chip that has
AVX2 and FMA, which looked like an obvious oversight.  It is not:

    -march=native -ffp-contract=off   58.1 s    0.0%   byte-identical
    -march=native                     57.3 s    1.5%   field total moves 2.3e-5
    -march=native -fassociative-math  55.3 s    5.0%   field total moves 7e-7..6.7e-5
    -march=native -ffast-math         54.4 s    6.5%   field total moves 3.6e-5

Widening the registers alone buys nothing, because the hot code is scalar libm
calls and pointer chasing through the stand and PFT lists, not loops a compiler
can vectorise.  Everything beyond that comes from reassociation, and costs bit
reproducibility.

How much it costs is worth stating precisely, because the headline number is
alarming and the reality is not.  Under `-ffast-math`, 8% of the GCGP values
change and the largest single relative change is 106% -- but the median change
is zero, the 99th percentile is 4.7e-4, and the field total moves by 3.6 parts
in 100,000 over three years.  The 106% is a value that passes through zero, so
its relative error is meaningless.  The model amplifies rounding differences
cell by cell while conserving the aggregate.

`-ffast-math` is nonetheless the wrong choice here, and not because of the
arithmetic.  It implies `-ffinite-math-only`, which permits GCC to fold
`isnan()` to false, and LPJmL has 64 `isnan`/`isinf` sites -- among them the
NetCDF climate readers, where they validate the input.  Trading away detection
of corrupt forcing data to save 1.5% over `-fassociative-math` is a bad bargain.

`-march=native -fassociative-math -fno-signed-zeros -fno-trapping-math
-freciprocal-math` keeps NaN semantics, gets 5.0%, and introduces no non-finite
value that was not already there.  It is left off by default: bit
reproducibility is what made the restart bug above findable, and 5% is not
worth losing it without the owner of the run deciding so.  To turn it on, add
those flags to `OPTFLAGS` in `Makefile.inc` after `./configure.sh`.

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
