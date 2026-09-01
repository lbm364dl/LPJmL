# Two outputs change under the cost-based split, and it is not understood

Tracked as issue #37. Found 2026-09-01 by comparing two runs of the same
1750–2023 global configuration that differ only in the decomposition and the
build; this file carries the detail and the reproduction recipe, the issue
carries the status.

**Nothing is broken downstream today** — neither field is consumed by WHEP. The
reason to chase it is that the mechanism, if it is what it looks like, is not
specific to these two fields.

## The observation

| | v1 (2026-08-27) | v2 (2026-09-01) |
|---|---|---|
| decomposition | default, equal cells per task | **cost-based** |
| build | plain | **PGO + LTO** |
| ranks | 30, unpinned | **32, pinned, huge pages** |
| code | `7c7ff9ed` | `bedc436e` |

Same forcing, same spinup, `min_cropfrac` unset in both. Field by field at seven
years spanning 1750–2023:

    fields compared : 169
      bit-identical : 166
      differing     : 3

`cft_airrig_month` differs by design (`a819b5d0`). The other two do not.

**`soilc_1m`** — exactly 60 of 58,795 cells differ, scattered across latitudes,
and the values look *swapped in magnitude*:

    [45,598]   v2 11.826    v1  0.86986
    [49,415]   v2  1.3702   v1 17.614
    [55,666]   v2 16.147    v1  1.1922
    [62,641]   v2  0.50087  v1  9.2384

Every other cell bit-identical. Global sum moves 0.01–0.03%.

**`pft_agtop_litterc`** — 19–29% of entries differ, dominated by near-zero
values (v2 holds 4.4e-11 where v1 holds exactly 0). v1 also carries small
negative above-ground litter carbon, minimum −0.0162, which is nonphysical;
v2's minimum is 0. Global sum moves −0.7% to +2.5% by year.

## Why the split is the suspect

`cc2e8ed1 fix(netcdf): stop assuming every task holds the same number of cells`
establishes that the assumption existed and was fixed in at least one place.
Under the default split every task holds the same count and it is harmless;
under the cost-based split it is not. Values landing in the wrong cells with
everything else exact is what a surviving instance would produce.

This is a hypothesis. It has not been tested, and **which run is correct for
these two fields is not established** — v2 having no negative litter carbon is
weak evidence in its favour and nothing more.

## Why it matters beyond these two fields

The 166 bit-identical fields are strong evidence the model state is unaffected:
`soilc`, `litfallc`, `vegc` and `npp` all match exactly. So this lives in the
output path, not the physics. But an output path that can misplace a value for
two fields can presumably do it for others, and the only reason it was visible
is that a run existed with the old split to compare against. Every future run
will use the cost-based split, and then there is nothing to compare with.

Both affected fields are aggregated or derived rather than direct state, which
may be the common factor.

## Cheapest first step

`soilc_1m`: 60 cells, exact swaps rather than perturbations, so following one
cell through `writepft` / `write_pft_float_netcdf` under an unequal-cell
decomposition should show it. Reproduce with any short run of the same
configuration with and without `bench/enable_balance.py --cost ...`, comparing
those two fields — a 2-year restart-based run on the full grid is about four
minutes.

Comparison tooling used: `compare_runs.py` in the v2 run's `configurations/`.
