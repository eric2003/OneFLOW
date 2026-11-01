#!/usr/bin/env python3
from __future__ import annotations
import argparse, csv
from pathlib import Path

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    xs, ys = [], []
    with args.csv.open(newline="", encoding="ascii") as fh:
        for row in csv.DictReader(fh):
            xs.append(float(row["chi"]))
            ys.append(float(row["fv1"]))
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(xs, ys, lw=2)
    ax.set_xlabel(r"$\chi$")
    ax.set_ylabel(r"$f_{v1}$")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=120)

if __name__ == "__main__":
    main()
