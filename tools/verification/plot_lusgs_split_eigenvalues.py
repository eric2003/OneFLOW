#!/usr/bin/env python3
"""Plot LUSGS split eigenvalues vs matrix sign (fixed flow state)."""
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
    s, l1, l2, l3 = [], [], [], []
    with args.csv.open(newline="", encoding="ascii") as fh:
        for row in csv.DictReader(fh):
            s.append(float(row["sign"]))
            l1.append(float(row["lmd1"]))
            l2.append(float(row["lmd2"]))
            l3.append(float(row["lmd3"]))
    fig, ax = plt.subplots(figsize=(6.5, 4.0))
    ax.plot(s, l1, 'o-', label=r"$\lambda_1$")
    ax.plot(s, l2, 's-', label=r"$\lambda_2$")
    ax.plot(s, l3, '^-', label=r"$\lambda_3$")
    ax.set_xlabel(r"matrix sign $s$")
    ax.set_ylabel(r"split $\lambda$")
    ax.set_title(r"LUSGS face eigenvalue split (vn=1, c=2, $|\lambda|_f=0.4$)")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best")
    fig.tight_layout()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=140)

if __name__ == "__main__":
    main()
