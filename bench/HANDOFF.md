# Running LPJmL fast on this machine

Written at the end of the performance work of 2026-08-28..31, for whoever picks
it up next.  `bench/README.md` has the reasoning and the measurements; this file
has the recipe and the state.

## The recipe

Three things, applied together.  All are byte-identical; `min_cropfrac` is a
fourth that is not, and is discussed separately below.

```sh
# 1. split the grid by measured cost instead of by cell count
bench/enable_balance.py <config.json> --cost /home/usuario/lpjml_perf_runs/cost_pcore.bin

# 2. build with a profile
bench/build_pgo.sh <a representative config.json>

# 3. run on all 32 hardware threads, pinned, with huge pages
mpirun -x GLIBC_TUNABLES=glibc.malloc.hugetlb=1 -np 32 \
       --bind-to hwthread --map-by hwthread ...
```

Also set `"write_cellcost_filename"` in the config.  The cost file accumulates
over every year simulated, so a full run with it set produces a
period-averaged file for the next run, which is worth about 2% more than the
1901-measured one.

## What each part is worth

Measured on `config_recover.json` -- global grid, river routing, 1901-2023
restart-based -- at one rank where repeatability is 0.1 s, and at 32 ranks in
both orders:

| | |
|---|---|
| cost-based split, 32 ranks pinned, against the default 24 unpinned | the bulk of it |
| `-lto` and a profile | 4.2% at one rank, 3.2% at 32 |
| `GLIBC_TUNABLES=glibc.malloc.hugetlb=1` | 0.9% |
| compute work merged into the binary | 12.7% at one rank |
| **all together, end to end, five years, both orders** | **520.7 s -> 288.9 s, 1.80x, byte-identical** |

Roughly a third of that is the binary and two thirds the decomposition.

## min_cropfrac: depends on how the run starts

Folds crop areas below a threshold into the largest crop of the same irrigation
regime.  The land-use input gives every crop a country reports a share in every
cropland cell of that country, so a fifth of all crop entries are under one
hectare and the smallest is a denormal float (reported upstream as
eduaguilera/whep#985).

  * **Started from scratch with a spinup** -- `config_rerun_recommended.json`,
    1450 to 2023 with output from 1750 -- there is no transient and it is worth
    **19%** on that shape.  First output year identical, then a steady -1.7% in
    nitrogen losses.
  * **Started from a restart written without it** -- `config_recover.json` --
    the stands it removes dump their nitrogen and losses run 54% high for six
    years before settling.  Leave it off, or discard those years.

Default is 0.  Set it per run: `bench/enable_balance.py <config> --min-cropfrac 1e-4`.

## Two switches deliberately left off

  * `./configure.sh -nosafe` drops 78 internal consistency checks.  4%, and
    byte-identical while none of them would have fired -- but a run that goes
    wrong then does so silently.
  * `-march=native -fassociative-math -fno-signed-zeros -fno-trapping-math
    -freciprocal-math` is a further 4% and gives up bit reproducibility.  Its
    perturbation is 260 times smaller than the model's own lambda-solver
    spread, so the arithmetic case is fine; the reason to decline is that
    byte-identical output is a working regression test and it is how the
    `fallow` bug was found.

## Measuring anything here

  * **Pin the ranks.**  Pinned 32-rank runs repeat to 0.15%; unpinned ones
    spread 12%.  Every early measurement in this work was fighting that, and two
    conclusions had to be reversed because of it.
  * **One rank is the honest place** to measure per-cell work: 0.1 s
    repeatability, and the 111 W package limit does not bite.
  * **`make clean` is not enough between experiments.**  `Makefile.inc` survives
    it, so profiler or other flags carry into the next build.  Always
    `make clean && ./configure.sh ...`.
  * **Profile more than one grid block.**  In 5000-5999 `littersom` is 6% of the
    run and the output path is 10%; in 36000-36999 it is 25% and 2%.
  * **Check nitrogen, not just carbon and water.**  `min_cropfrac` moved carbon
    pools not at all and nitrogen losses by half.  The columns are in
    `output/run.out`.

## Artefacts on disk

    /home/usuario/lpjml_perf_runs/cost_pcore.bin     per-cell cost, measured on
                                                     8 P-cores so core speed is
                                                     not baked in
    /home/usuario/lpjml_perf_runs/weights_24.csv     task weights, 24 ranks only;
                                                     unnecessary at 32
    /home/usuario/lpjml_perf_runs/cost_period_REJECTED_do_not_use.bin

## State of the fork

Everything from this work is merged to `master` (PRs #8 through #33).  Three
branches carry commits that are not:

  * `fix/litterfall-outputs-and-decay` -- **two output fixes by the repository
    owner, dated 2026-08-31, not in master.**  Not part of this work; worth a
    look.
  * `perf/multichunk` and `perf/lpjml-speedup` -- an abandoned multi-chunk
    decomposition and a multi-band cost-file format.  `bench/README.md` records
    why chunking was not pursued.  Nothing depends on either.
  * `rescue/dirty-tree-20260818` -- predates this work.
