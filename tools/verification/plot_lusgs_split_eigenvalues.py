#!/usr/bin/env python3
"""Plot LUSGS split eigenvalues vs matrix sign."""
from __future__ import annotations
import argparse, csv, sys
from pathlib import Path

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    tools = Path(__file__).resolve().parent
    if str(tools) not in sys.path:
        sys.path.insert(0, str(tools))
    from plot_style import apply_theory_style, save_figure
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
    fig, ax = plt.subplots(figsize=(7.0, 4.2))
    ax.plot(s, l1, "o-", lw=2, label="lambda1 (vn branch)")
    ax.plot(s, l2, "s-", lw=2, label="lambda2 (vn+c branch)")
    ax.plot(s, l3, "^-", lw=2, label="lambda3 (vn-c branch)")
    apply_theory_style(
        ax,
        xlabel="matrix sign s",
        ylabel="split lambda",
        title="LUSGS face eigenvalue split vs matrix sign",
        subtitle="ref: vn=1, c=2, |lambda|_f=0.4;  lambda=0.5*(phys + s*|lambda|_f)",
        legend=True,
    )
    save_figure(fig, args.out)
    print(f"wrote {args.out}")

if __name__ == "__main__":
    main()
