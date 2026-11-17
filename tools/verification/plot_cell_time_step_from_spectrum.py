#!/usr/bin/env python3
"""Local cell time step vs spectral radius."""
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
            xs.append(float(row["spectral_radius"]))
            ys.append(float(row["dt"]))
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    ax.plot(xs, ys, lw=2.0, label=r"$\Delta t=\mathrm{CFL}\,V/\lambda$")
    if apply_theory_style:
        apply_theory_style(
            ax,
            xlabel=r"spectral radius $\lambda$",
            ylabel=r"$\Delta t$",
            title=r"Cell time step from spectrum",
            subtitle=r"ref $\mathrm{CFL}=1$, $V=2$: $\Delta t=\mathrm{CFL}\,V/\lambda$",
        )
    else:
        ax.set_xlabel(r"lambda")
        ax.set_ylabel(r"dt")
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
