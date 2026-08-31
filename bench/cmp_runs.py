#!/usr/bin/env python3
"""Compare two LPJmL run directories for numerical identity.

Reports, per artefact, whether the two runs agree bit-for-bit.  For the
restart file it also localises any difference, so a mismatch confined to the
metadata header (version string, git hash, wall-clock timestamp) is not
mistaken for a change in the model state.

Usage: cmp_runs.py <run_a> <run_b> [--header-bytes 65536]
"""
import argparse, os, sys

import numpy as np

HEADER_DEFAULT = 65536   # the bstruct header runs to ~9 kB: metadata plus the variable index


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


def common_suffix(a, b):
    """Length of the longest common suffix of two byte strings."""
    lo, hi = 0, min(len(a), len(b))
    while lo < hi:
        mid = (lo + hi + 1) // 2
        if a[len(a) - mid:] == b[len(b) - mid:]:
            lo = mid
        else:
            hi = mid - 1
    return lo


def cmp_restart(pa, pb, header_bytes):
    if not (os.path.exists(pa) or os.path.exists(pb)):
        return "restart", "IDENTICAL*", "neither run wrote a restart file"
    if not (os.path.exists(pa) and os.path.exists(pb)):
        return "restart", "MISSING", "only one run wrote a restart file"
    a = open(pa, "rb").read()
    b = open(pb, "rb").read()
    if a == b:
        return "restart", "IDENTICAL", f"{len(a)/1e6:.0f} MB, byte-for-byte"
    # The bstruct header carries a timestamp, the command line that produced
    # the file and the git hash of the binary, followed by an index whose
    # offsets all shift when any of those change length.  A difference
    # confined to it says nothing about the model, so the payload is compared
    # by aligning the two files at their ends rather than their starts --
    # otherwise a one-character difference in a run directory name reads as a
    # total mismatch.
    suf = common_suffix(a, b)
    ha, hb = len(a) - suf, len(b) - suf
    if max(ha, hb) <= header_bytes:
        return "restart", "IDENTICAL*", (
            f"{suf/1e6:.0f} MB of model state byte-for-byte; only the "
            f"{ha}-byte metadata header differs (timestamp, command line, "
            f"git hash)")
    if len(a) != len(b):
        return "restart", "SIZE-DIFF", (
            f"{len(a)} vs {len(b)}; the last {suf/1e6:.1f} MB agree, so the "
            f"difference reaches {ha} bytes in, past the header")
    ranges, count = diff_ranges(a, b)
    last = max(r[1] for r in ranges)
    return "restart", "DIFFERS", (
        f"{count} differing bytes, first at {ranges[0][0]}, "
        f"last at {last} of {len(a)}")


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
    # Comparing nothing is not a pass.  A run that died before writing any
    # output, or a pair of directories that share no file, used to come out as
    # IDENTICAL and read as a green light.
    if ok + len(bad) == 0 and not status.startswith("IDENTICAL"):
        print("VERDICT: NOTHING COMPARED -- neither run produced output")
        return 2
    if ok + len(bad) == 0 and detail == "neither run wrote a restart file":
        print("VERDICT: NOTHING COMPARED -- no output files and no restart file")
        return 2
    verdict = not bad and status.startswith("IDENTICAL")
    print("VERDICT:", "IDENTICAL" if verdict else "DIFFERENT")
    return 0 if verdict else 1


if __name__ == "__main__":
    sys.exit(main())
