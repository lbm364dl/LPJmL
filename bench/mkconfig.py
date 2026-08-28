#!/usr/bin/env python3
"""Derive a benchmark LPJmL config from the production WHEP config.

The production config is the one that produced the 11h29m global run
(30 ranks, 300 spinup years + 1901-2023).  Benchmarks keep every physical
switch identical and only shrink the amount of work, so that whatever we
measure is the same code path the real run takes.
"""
import json, os, sys, argparse

# The 71-output production config as it stands on the current master.  The
# original in lpjml_rerun_2026-08 predates two merged PRs that added 13 output
# variables, so LPJmL rejects its 384-entry outputvar table against NOUT=397.
PROD = "/home/usuario/lpjml_perf_runs/config_prod_397.json"
OLD_PREFIX = "/home/usuario/LPJmL-611/global_1901-2023_spinup_300_611"


def retarget(obj, old, new):
    """Rewrite every absolute path that points into the old run directory."""
    if isinstance(obj, str):
        return obj.replace(old, new)
    if isinstance(obj, list):
        return [retarget(v, old, new) for v in obj]
    if isinstance(obj, dict):
        return {k: retarget(v, old, new) for k, v in obj.items()}
    return obj


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True, help="run directory to create")
    ap.add_argument("--base", default=PROD)
    ap.add_argument("--old-prefix", dest="old_prefix", default=OLD_PREFIX,
                    help="run directory the base config points at, rewritten to --out")
    ap.add_argument("--startgrid", default=None)
    ap.add_argument("--endgrid", default=None)
    ap.add_argument("--nspinup", type=int, default=None)
    ap.add_argument("--firstyear", type=int, default=None)
    ap.add_argument("--lastyear", type=int, default=None)
    ap.add_argument("--river-routing", dest="rr", default=None,
                    choices=["true", "false"])
    ap.add_argument("--fmt", default=None, help="default_fmt override, e.g. raw")
    ap.add_argument("--set", action="append", default=[],
                    help="extra top-level override, key=json_value")
    args = ap.parse_args()

    cfg = json.load(open(args.base))
    run = os.path.abspath(args.out)
    cfg = retarget(cfg, args.old_prefix, run)

    if args.startgrid is not None:
        cfg["startgrid"] = int(args.startgrid) if args.startgrid != "all" else "all"
    if args.endgrid is not None:
        cfg["endgrid"] = int(args.endgrid) if args.endgrid != "all" else "all"
    if args.nspinup is not None:
        cfg["nspinup"] = args.nspinup
    if args.firstyear is not None:
        cfg["firstyear"] = args.firstyear
        cfg["outputyear"] = args.firstyear
    if args.lastyear is not None:
        cfg["lastyear"] = args.lastyear
    if args.rr is not None:
        cfg["river_routing"] = (args.rr == "true")
    if args.fmt is not None:
        cfg["default_fmt"] = args.fmt
        for o in cfg["output"]:
            f = o["file"]
            if f.get("fmt") not in ("txt",):
                f.pop("fmt", None)
    for kv in args.set:
        # value is JSON so that strings, numbers, booleans and null all work
        k, v = kv.split("=", 1)
        cfg[k] = json.loads(v)

    os.makedirs(os.path.join(run, "configurations"), exist_ok=True)
    os.makedirs(os.path.join(run, "output", "scenario_1"), exist_ok=True)
    os.makedirs(os.path.join(run, "restart", "scenario_1"), exist_ok=True)
    dest = os.path.join(run, "configurations", "config_scenario_1.json")
    with open(dest, "w") as fp:
        json.dump(cfg, fp, indent=1)
    print(dest)


if __name__ == "__main__":
    main()
