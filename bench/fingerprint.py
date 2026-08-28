#!/usr/bin/env python3
"""Correctness fingerprint of an LPJmL run directory.

Two runs that are numerically identical must produce the same fingerprint.
Covers:
  * every NetCDF output variable's raw data bytes (metadata excluded, because
    LPJmL stamps creation time and command line into the file);
  * every raw/clm/txt output byte-for-byte;
  * the restart file from the end of its header onward (the header carries a
    wall-clock timestamp, the rest is the complete model state).

Usage: fingerprint.py <run_dir> [--json out.json]
"""
import argparse, hashlib, json, os, sys

import numpy as np


def sha(b):
    return hashlib.sha256(b).hexdigest()[:32]


def restart_fingerprint(path):
    """Hash the restart file past its header.

    The bstruct header ends with the key of the cell index array, "grid".
    Everything before it is metadata (version string, git hash, and a
    history string containing the run's wall-clock timestamp).
    """
    with open(path, "rb") as fp:
        head = fp.read(65536)
        i = head.rfind(b"grid")
        if i < 0:
            return {"error": "grid marker not found", "size": os.path.getsize(path)}
        fp.seek(i)
        h = hashlib.sha256()
        while True:
            chunk = fp.read(1 << 22)
            if not chunk:
                break
            h.update(chunk)
    return {"offset": i, "size": os.path.getsize(path), "sha": h.hexdigest()[:32]}


def netcdf_fingerprint(path):
    from netCDF4 import Dataset
    out = {}
    with Dataset(path) as ds:
        for name, var in ds.variables.items():
            if name in ("lon", "lat", "time", "longitude", "latitude"):
                continue
            data = np.ma.filled(var[:], np.nan)
            arr = np.ascontiguousarray(data)
            out[name] = {"shape": list(arr.shape), "sha": sha(arr.tobytes())}
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run_dir")
    ap.add_argument("--json", default=None)
    args = ap.parse_args()

    fp = {}
    outdir = os.path.join(args.run_dir, "output", "scenario_1")
    for name in sorted(os.listdir(outdir)) if os.path.isdir(outdir) else []:
        path = os.path.join(outdir, name)
        if not os.path.isfile(path):
            continue
        if name.endswith(".json") or name.endswith(".log"):
            continue          # metadata sidecars carry timestamps
        try:
            if name.endswith(".nc") or name.endswith(".nc4"):
                fp["nc:" + name] = netcdf_fingerprint(path)
            elif name.endswith(".csv") or name.endswith(".txt"):
                fp["txt:" + name] = sha(open(path, "rb").read())
            else:
                fp["bin:" + name] = sha(open(path, "rb").read())
        except Exception as exc:                       # noqa: BLE001
            fp["err:" + name] = repr(exc)

    rst = os.path.join(args.run_dir, "restart", "scenario_1", "restart.lpj")
    if os.path.exists(rst):
        fp["restart"] = restart_fingerprint(rst)

    text = json.dumps(fp, indent=1, sort_keys=True)
    if args.json:
        open(args.json, "w").write(text)
    print(text)


if __name__ == "__main__":
    main()
