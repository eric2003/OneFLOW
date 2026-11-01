#!/usr/bin/env python3
"""Plot SST F2(arg2)=tanh(arg2^2) theory curve (read-only CSV under repo)."""
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
            xs.append(float(row["arg2"]))
            ys.append(float(row["F2"]))
    fig, ax = plt.subplots(figsize=(6.5, 4.0))
    ax.plot(xs, ys, lw=2.0, color="C1", label=r"$F_2=\tanh(\mathrm{arg}_2^2)$")
    ax.axhline(1.0, color="0.5", ls="--", lw=1)
    ax.set_xlabel(r"$\mathrm{arg}_2$")
    ax.set_ylabel(r"$F_2$")
    ax.set_title(r"SST $F_2(\mathrm{arg}_2)$ intermediate theory curve")
    ax.set_xlim(0, 4)
    ax.set_ylim(0, 1.05)
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best")
    fig.tight_layout()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=140)

if __name__ == "__main__":
    main()
