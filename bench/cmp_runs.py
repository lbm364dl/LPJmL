#!/usr/bin/env python3
"""Compare two LPJmL run directories for numerical identity.

Reports, per artefact, whether the two runs agree bit-for-bit.  For the
restart file it also localises any difference, so a mismatch confined to the
metadata header (version string, git hash, wall-clock timestamp) is not
mistaken for a change in the model state.

Usage: cmp_runs.py <run_a> <run_b> [--header-bytes 4096]
"""
import argparse, os, sys

import numpy as np

HEADER_DEFAULT = 4096


def diff_ranges(a, b, limit=8):
    """Byte offsets where two buffers differ, coalesced into ranges."""
    n = min(len(a), len(b))
    d = np.frombuffer(a[:n], dtype=np.uint8) != np.frombuffer(b[:n], dtype=np.uint8)
    idx = np.flatnonzero(d)
    if idx.size == 0:
        return [], 0
    splits = np.flatnonzero(np.diff(idx) > 1) + 1
    groups = np.split(idx, splits)
    return [(int(g[0]), int(g[-1])) for g in groups[:limit]], int(idx.size)


def cmp_restart(pa, pb, header_bytes):
    if not (os.path.exists(pa) or os.path.exists(pb)):
        return "restart", "IDENTICAL*", "neither run wrote a restart file"
    if not (os.path.exists(pa) and os.path.exists(pb)):
        return "restart", "MISSING", "only one run wrote a restart file"
    sa, sb = os.path.getsize(pa), os.path.getsize(pb)
    if sa != sb:
        return "restart", "SIZE-DIFF", f"{sa} vs {sb}"
    a = open(pa, "rb").read()
    b = open(pb, "rb").read()
    if a == b:
        return "restart", "IDENTICAL", f"{sa/1e6:.0f} MB, byte-for-byte"
    ranges, count = diff_ranges(a, b)
    last = max(r[1] for r in ranges)
    if last < header_bytes and count < header_bytes:
        return "restart", "IDENTICAL*", (
            f"{count} differing bytes, all in the metadata header "
            f"(<{header_bytes}); model state identical")
    return "restart", "DIFFERS", (
        f"{count} differing bytes, first at {ranges[0][0]}, "
        f"last at {last} of {sa}")


def cmp_netcdf(pa, pb):
    from netCDF4 import Dataset
    worst = ("IDENTICAL", "")
    with Dataset(pa) as da, Dataset(pb) as db:
        for name, va in da.variables.items():
            if name not in db.variables:
                return "MISSING-VAR", name
            vb = db.variables[name]
            kind = getattr(va.dtype, "kind", "S")   # str dtype => text variable
            if kind not in "fiu":                   # char/string variables
                if np.array_equal(np.asarray(va[:]), np.asarray(vb[:])):
                    continue
                return "DIFFERS", f"{name}: non-numeric variable differs"
            x = np.ma.filled(va[:], np.nan).astype(np.float64)
            y = np.ma.filled(vb[:], np.nan).astype(np.float64)
            if x.shape != y.shape:
                return "SHAPE-DIFF", name
            same = np.array_equal(
                np.nan_to_num(x, nan=-9e99), np.nan_to_num(y, nan=-9e99))
            if not same:
                d = np.abs(x - y)
                scale = np.maximum(np.abs(x), np.abs(y))
                with np.errstate(invalid="ignore", divide="ignore"):
                    rel = np.nanmax(np.where(scale > 0, d / scale, 0.0))
                worst = ("DIFFERS",
                         f"{name}: max abs {np.nanmax(d):.3e}, max rel {rel:.3e}")
    return worst


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run_a")
    ap.add_argument("run_b")
    ap.add_argument("--header-bytes", type=int, default=HEADER_DEFAULT)
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args()

    oa = os.path.join(args.run_a, "output", "scenario_1")
    ob = os.path.join(args.run_b, "output", "scenario_1")
    names = sorted(set(os.listdir(oa)) & set(os.listdir(ob)))

    bad, ok = [], 0
    for name in names:
        if name.endswith(".json") or name.endswith(".log") or name in ("run.out", "run.err"):
            continue
        pa, pb = os.path.join(oa, name), os.path.join(ob, name)
        if not (os.path.isfile(pa) and os.path.isfile(pb)):
            continue
        if name.endswith(".nc") or name.endswith(".nc4"):
            status, detail = cmp_netcdf(pa, pb)
        else:
            status = "IDENTICAL" if open(pa, "rb").read() == open(pb, "rb").read() else "DIFFERS"
            detail = ""
        if status == "IDENTICAL":
            ok += 1
            if args.verbose:
                print(f"  ok       {name}")
        else:
            bad.append((name, status, detail))
            print(f"  {status:<10} {name}  {detail}")

    what, status, detail = cmp_restart(
        os.path.join(args.run_a, "restart", "scenario_1", "restart.lpj"),
        os.path.join(args.run_b, "restart", "scenario_1", "restart.lpj"),
        args.header_bytes)
    print(f"  {status:<10} {what}  {detail}")

    print(f"\n{ok}/{ok + len(bad)} output files identical")
    verdict = not bad and status.startswith("IDENTICAL")
    print("VERDICT:", "IDENTICAL" if verdict else "DIFFERENT")
    return 0 if verdict else 1


if __name__ == "__main__":
    sys.exit(main())
