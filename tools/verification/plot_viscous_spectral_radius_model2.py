#!/usr/bin/env python3
"""Viscous face spectral radius model2 vs total viscosity."""
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
    try:
        from plot_style import apply_theory_style, save_figure
    except Exception:
        apply_theory_style = None
        save_figure = None
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    xs, ys = [], []
    with args.csv.open(newline="", encoding="ascii") as fh:
        for row in csv.DictReader(fh):
            xs.append(float(row["total_viscosity"]))
            ys.append(float(row["lambda_vis"]))
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    ax.plot(xs, ys, lw=2.0, label=r"$\lambda_{\mathrm{vis}}$ (model 2)")
    if apply_theory_style:
        apply_theory_style(
            ax,
            xlabel=r"$\mu+\mu_t$",
            ylabel=r"$\lambda_{\mathrm{vis}}$",
            title=r"Viscous face spectral radius (model 2)",
            subtitle=r"$\lambda=0.5\cdot(2\mu_{tot})/(\rho d Re)\cdot A$; ref $d=0.1$, $Re=100$",
        )
    else:
        ax.set_xlabel(r"total viscosity")
        ax.set_ylabel(r"lambda_vis")
        ax.grid(True, alpha=0.35)
        ax.legend(loc="best")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    if save_figure:
        save_figure(fig, args.out)
    else:
        fig.tight_layout()
        fig.savefig(args.out, dpi=140, bbox_inches="tight", facecolor="white")
    print(f"wrote {args.out}")

if __name__ == "__main__":
    main()
