#!/usr/bin/env python3
"""Plot NS inviscid face spectral radius vs |vn| (A=2,c=1)."""
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
            xs.append(float(row["abs_vn"]))
            ys.append(float(row["lambda_inv"]))
    fig, ax = plt.subplots(figsize=(6.5, 4.0))
    ax.plot(xs, ys, lw=2.0, label=r"$\lambda=\frac{A}{2}(|v_n|+c)$")
    ax.set_xlabel(r"$|v_n|$")
    ax.set_ylabel(r"$\lambda_{\mathrm{inv}}$")
    ax.set_title("Mean-flow inviscid face spectral radius (A=2, c=1)")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best")
    fig.tight_layout()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=140)

if __name__ == "__main__":
    main()
