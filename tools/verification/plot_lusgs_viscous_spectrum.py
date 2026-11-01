#!/usr/bin/env python3
"""Plot LUSGS viscous spectrum contribution vs dq (fixed visrad)."""
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
            xs.append(float(row["dq"]))
            ys.append(float(row["contrib"]))
    fig, ax = plt.subplots(figsize=(6.5, 4.0))
    ax.plot(xs, ys, lw=2.0, color="C6",
            label=r"$\Delta r=-\lambda_{\mathrm{vis}}\Delta q$ ($\lambda=2$)")
    ax.axhline(0.0, color="0.5", ls="--", lw=1)
    ax.axvline(0.0, color="0.5", ls=":", lw=1)
    ax.set_xlabel(r"$\Delta q$")
    ax.set_ylabel(r"$\Delta r$")
    ax.set_title("LUSGS viscous spectral residual contribution")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best")
    fig.tight_layout()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=140)

if __name__ == "__main__":
    main()
