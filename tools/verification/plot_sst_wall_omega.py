#!/usr/bin/env python3
"""Plot Menter wall omega vs wall distance (theory)."""
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
            xs.append(float(row["d"]))
            ys.append(float(row["omega_w"]))
    fig, ax = plt.subplots(figsize=(6.5, 4.0))
    ax.plot(xs, ys, lw=2.0, color="C4", label=r"$\omega_w(d)$")
    ax.set_xlabel(r"$d$ (wall distance)")
    ax.set_ylabel(r"$\omega_w$")
    ax.set_title(r"Menter wall $\omega_w=60\mu/(\rho\beta_1 d^2 Re)$ (ref constants)")
    ax.set_xlim(0.05, 2.0)
    ax.set_yscale("log")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(loc="best")
    fig.tight_layout()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=140)

if __name__ == "__main__":
    main()
