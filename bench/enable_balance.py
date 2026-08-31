#!/usr/bin/env python3
"""Turn on the cost-based decomposition in an LPJmL configuration.

    enable_balance.py <config.json> --cost cost_pcore.bin [--weights w.csv|CSV]
    enable_balance.py <config.json> --min-cropfrac 1e-4
    enable_balance.py <config.json> --off

A configuration written by lpjmlkit carries its own copy of the parameters
rather than reading par/lpjparam.cjson, so a parameter changed in the model
tree does not reach a configuration already generated.  --min-cropfrac writes
it into the copy.

The cost file is produced by any run with "write_cellcost_filename" set, and
should be measured on cores of one kind -- eight P-cores, say -- because the
file records what each cell cost on the rank that owned it, and a run spread
over a hybrid chip charges E-core cells for the core rather than the cell.

The weights say how fast each task is, and only mean anything if the ranks are
pinned:

    mpirun -np 24 --bind-to core --map-by core ...

Derive them with cellcost.py retune from a pinned run that wrote its own cost
file, and re-derive them if the rank count, the machine or the binary changes.
A stale cost file or stale weights cost performance and never accuracy: they can
only produce a worse split, never a wrong answer.
"""
import argparse
import json
import os
import sys


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("config")
    ap.add_argument("--cost", help="per-cell cost file to split by")
    ap.add_argument("--weights", default=None,
                    help="comma-separated per-task speeds, or a file holding them")
    ap.add_argument("--off", action="store_true", help="remove both keys again")
    ap.add_argument("--min-cropfrac", dest="min_cropfrac", type=float, default=None,
                    help="fold crop areas below this fraction of the cell into "
                         "the largest crop of the same irrigation regime")
    args = ap.parse_args()

    with open(args.config) as f:
        cfg = json.load(f)

    if args.off:
        for k in ("cellcost_filename", "task_weights"):
            cfg.pop(k, None)
        cfg.get("param", {}).pop("min_cropfrac", None)
        print("removed cellcost_filename, task_weights and min_cropfrac")
    elif args.min_cropfrac is not None and not args.cost:
        cfg.setdefault("param", {})["min_cropfrac"] = args.min_cropfrac
        print(f"param.min_cropfrac = {args.min_cropfrac}")
    else:
        if not args.cost:
            sys.exit("--cost or --min-cropfrac is required unless --off is given")
        cost = os.path.abspath(args.cost)
        if not os.path.exists(cost):
            sys.exit(f"no cost file at {cost}")
        cfg["cellcost_filename"] = cost
        print(f"cellcost_filename = {cost}")
        if args.weights:
            text = args.weights
            if os.path.exists(text):
                text = open(text).read()
            w = [float(x) for x in text.replace("[", "").replace("]", "")
                 .replace("\n", ",").split(",") if x.strip()]
            cfg["task_weights"] = w
            print(f"task_weights = {len(w)} values, "
                  f"{min(w):.3f}..{max(w):.3f}")
            print("  these only apply if the ranks are pinned: "
                  "mpirun --bind-to core --map-by core")
        else:
            cfg.pop("task_weights", None)
        if args.min_cropfrac is not None:
            cfg.setdefault("param", {})["min_cropfrac"] = args.min_cropfrac
            print(f"param.min_cropfrac = {args.min_cropfrac}")

    with open(args.config, "w") as f:
        json.dump(cfg, f, indent=1)
    print(f"wrote {args.config}")


if __name__ == "__main__":
    main()
