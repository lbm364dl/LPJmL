#!/usr/bin/env python3
"""Build, inspect and evaluate LPJmL cell cost files.

The model splits the grid into equal-sized contiguous blocks, one per task.
The grid runs north to south and cell cost spans more than an order of
magnitude, so equal counts are not equal work; with river routing the tasks
meet at a barrier every simulated day and the run goes at the pace of the
heaviest task.  `cellcost_filename` in the config makes the split follow cost
instead.

The authoritative cost file comes from a calibration run
(`write_cellcost_filename`), which times every cell.  This tool can also build
a proxy from the land-use input when no calibration run is available, and can
report what imbalance any cost file implies.

    cellcost.py proxy   --out cost.bin [--base 1 --per-stand 0.3]
    cellcost.py dump    cost.bin
    cellcost.py report  cost.bin --ntask 30 [--weights w0,w1,...]

File format (little-endian, as written by writecellcost() in
src/lpj/divide_cells.c): "LPJCOST", int32 version=1, int32 firstcell,
int32 ncell, then ncell float64 costs in cell order.
"""
import argparse, struct, sys

import numpy as np

HEADER = b"LPJCOST"
VERSION = 1
INPUTS = "/home/usuario/WHEP/LPJmL_inputs/whep/lpjml_inputs"
GRID = f"{INPUTS}/gadm/grid_gadm_30arcmin.bin"
LANDUSE = f"{INPUTS}/landuse/cft_default_cft_aggregation_30min_1851-2023.nc"


def read_grid(path=GRID):
    raw = open(path, "rb").read()
    if raw[:7] != b"LPJGRID":
        raise SystemExit(f"{path} is not an LPJGRID file")
    ver, order, firstyear, nyear, firstcell, ncell, nbands = struct.unpack_from("<7i", raw, 7)
    scalar = struct.unpack_from("<f", raw, 7 + 32)[0]
    xy = np.frombuffer(raw, "<i2", count=ncell * 2, offset=7 + 44).astype(np.float64) * scalar
    return firstcell, xy[0::2], xy[1::2]


def write_cost(path, firstcell, cost):
    with open(path, "wb") as fp:
        fp.write(HEADER)
        fp.write(struct.pack("<3i", VERSION, firstcell, cost.size))
        fp.write(np.ascontiguousarray(cost, dtype="<f8").tobytes())
    print(f"wrote {cost.size} cell costs to {path}")


def read_cost(path):
    raw = open(path, "rb").read()
    if raw[:7] != HEADER:
        raise SystemExit(f"{path} is not an LPJCOST file")
    ver, firstcell, ncell = struct.unpack_from("<3i", raw, 7)
    if ver != VERSION:
        raise SystemExit(f"{path} has version {ver}, expected {VERSION}")
    cost = np.frombuffer(raw, "<f8", count=ncell, offset=7 + 12)
    return firstcell, cost


def divide_equal(n, ntask):
    """Exact port of the historical divide() in src/lpj/fscanconfig.c."""
    bounds, lo = [0], 0
    hi = n // ntask - 1 + (1 if n % ntask else 0)
    for i in range(1, ntask):
        bounds.append(hi + 1)
        lo = hi + 1
        hi = lo + n // ntask - 1
        if n % ntask > i:
            hi += 1
    bounds.append(n)
    return np.array(bounds)


def divide_cost(cost, ntask, weight=None):
    """Port of divide_cells() in src/lpj/divide_cells.c."""
    n = cost.size
    w = np.ones(ntask) if weight is None else np.asarray(weight, float)
    c = np.cumsum(np.maximum(cost, 1e-9))
    targets = c[-1] * np.cumsum(w)[:-1] / w.sum()
    b = np.searchsorted(c, targets) + 1
    b = np.concatenate([[0], b, [n]])
    for r in range(1, ntask):
        b[r] = max(b[r], b[r - 1] + 1)
        b[r] = min(b[r], n - (ntask - r))
    return b


def loads(cost, bounds, ntask, weight=None):
    w = np.ones(ntask) if weight is None else np.asarray(weight, float)
    return np.array([cost[bounds[i]:bounds[i + 1]].sum() for i in range(ntask)]) / w


def cmd_proxy(args):
    from netCDF4 import Dataset
    firstcell, lon, lat = read_grid()
    n = lon.size
    with Dataset(LANDUSE) as ds:
        flon = ds.variables["longitude"][:].data
        flat = ds.variables["latitude"][:].data
        t = ds.variables["time"][:].data
        ti = int(np.argmin(np.abs(t - args.year)))
        lu = np.ma.filled(ds.variables["landuse"][ti], 0.0)
    ilon = np.searchsorted(flon, lon - 1e-6)
    ilat = (len(flat) - 1 - np.searchsorted(flat[::-1], lat - 1e-6)
            if flat[0] > flat[-1] else np.searchsorted(flat, lat - 1e-6))
    ok = (ilat >= 0) & (ilat < len(flat)) & (ilon >= 0) & (ilon < len(flon))
    nstand = np.zeros(n)
    sel = lu[:, ilat[ok], ilon[ok]]
    sel[~np.isfinite(sel)] = 0.0
    nstand[ok] = (sel > 0).sum(axis=0)
    cost = args.base + args.per_stand * nstand
    print(f"proxy from land use {args.year}: stands per cell "
          f"min {nstand.min():.0f} mean {nstand.mean():.2f} max {nstand.max():.0f}")
    write_cost(args.out, firstcell, cost)


def cmd_dump(args):
    firstcell, cost = read_cost(args.file)
    print(f"firstcell={firstcell} ncell={cost.size} total={cost.sum():.1f}")
    q = np.percentile(cost, [0, 1, 25, 50, 75, 99, 100])
    print("cost percentiles  min/1/25/50/75/99/max: " + "  ".join(f"{v:.4g}" for v in q))
    print(f"heaviest cell / median cell = {cost.max() / max(np.median(cost), 1e-12):.1f}x")


def cmd_report(args):
    firstcell, cost = read_cost(args.file)
    n, ntask = cost.size, args.ntask
    weight = None
    if args.weights:
        weight = np.array([float(x) for x in args.weights.split(",")])
        if weight.size != ntask:
            raise SystemExit(f"{weight.size} weights for {ntask} tasks")
    for label, bounds in (("equal cell count (current)", divide_equal(n, ntask)),
                          ("by cost", divide_cost(cost, ntask, weight))):
        load = loads(cost, bounds, ntask, weight)
        sizes = np.diff(bounds)
        print(f"{label:28s} slowest/mean {load.max()/load.mean():6.3f}   "
              f"cells per task {sizes.min():5d}..{sizes.max():5d}   "
              f"idle {100*(1-load.mean()/load.max()):4.1f}%")
    base = loads(cost, divide_equal(n, ntask), ntask, weight).max()
    tuned = loads(cost, divide_cost(cost, ntask, weight), ntask, weight).max()
    print(f"\npredicted speedup from balancing alone: {base/tuned:.2f}x")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("proxy", help="build a cost file from the land-use input")
    p.add_argument("--out", required=True)
    p.add_argument("--year", type=int, default=1901)
    p.add_argument("--base", type=float, default=1.0,
                   help="cost of a cell with no crop stands")
    p.add_argument("--per-stand", type=float, default=0.3,
                   help="added cost per active crop stand")
    p.set_defaults(func=cmd_proxy)

    p = sub.add_parser("dump", help="summarise a cost file")
    p.add_argument("file")
    p.set_defaults(func=cmd_dump)

    p = sub.add_parser("report", help="imbalance implied by a cost file")
    p.add_argument("file")
    p.add_argument("--ntask", type=int, default=30)
    p.add_argument("--weights", default=None, help="comma-separated per-task speeds")
    p.set_defaults(func=cmd_report)

    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
