#!/usr/bin/env python3
"""Plot SST F1(arg1)=tanh(arg1^4) theory curve."""
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
            xs.append(float(row["arg1"]))
            ys.append(float(row["F1"]))
    fig, ax = plt.subplots(figsize=(6.5, 4.0))
    ax.plot(xs, ys, lw=2.0, color="C0", label=r"$F_1=\tanh(\mathrm{arg}_1^4)$")
    ax.axhline(1.0, color="0.5", ls="--", lw=1)
    ax.axvline(1.0, color="0.5", ls=":", lw=1)
    ax.set_xlabel(r"$\mathrm{arg}_1$")
    ax.set_ylabel(r"$F_1$")
    ax.set_title("Menter SST blending function $F_1$")
    ax.set_xlim(0, 4)
    ax.set_ylim(0, 1.05)
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best")
    fig.tight_layout()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=140)
    print(f"wrote {args.out}")

if __name__ == "__main__":
    main()
