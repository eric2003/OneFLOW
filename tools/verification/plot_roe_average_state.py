#!/usr/bin/env python3
"""Roe average density and velocity vs right density."""
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
    xs, dm, um = [], [], []
    with args.csv.open(newline="", encoding="ascii") as fh:
        for row in csv.DictReader(fh):
            xs.append(float(row["density_right"]))
            dm.append(float(row["density_mid"]))
            um.append(float(row["u_mid"]))
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    ax.plot(xs, dm, lw=2.0, label=r"$\rho_m=\sqrt{\rho_L\rho_R}$")
    ax.plot(xs, um, lw=2.0, ls="--", label=r"$u_m$ (ref $u_L=0$, $u_R=2$)")
    if apply_theory_style:
        apply_theory_style(
            ax,
            xlabel=r"$\rho_R$",
            ylabel=r"Roe mid state",
            title=r"Roe average state",
            subtitle=r"ref $\rho_L=1$",
        )
    else:
        ax.set_xlabel(r"rhoR")
        ax.set_ylabel(r"mid")
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
