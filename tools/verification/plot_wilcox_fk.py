#!/usr/bin/env python3
"""Plot Wilcox fk(xk) theory curve."""
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
            xs.append(float(row["xk"]))
            ys.append(float(row["fk"]))
    fig, ax = plt.subplots(figsize=(6.5, 4.0))
    ax.plot(xs, ys, lw=2.0, color="C3", label=r"$f_k(x_k)$")
    ax.axhline(1.0, color="0.5", ls="--", lw=1)
    ax.axvline(0.0, color="0.5", ls=":", lw=1)
    ax.set_xlabel(r"$x_k$")
    ax.set_ylabel(r"$f_k$")
    ax.set_title(r"Wilcox $f_k(x_k)$ (beta* modifier)")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best")
    fig.tight_layout()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=140)

if __name__ == "__main__":
    main()
