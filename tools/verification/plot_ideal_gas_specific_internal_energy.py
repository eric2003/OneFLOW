#!/usr/bin/env python3
"""Ideal-gas specific total internal energy vs kinetic energy."""
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
            xs.append(float(row["speed_squared"]))
            ys.append(float(row["specific_energy"]))
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    ax.plot(xs, ys, lw=2.0, label=r"$e=p/(\rho(\gamma-1))+\frac{1}{2}|u|^{2}$")
    if apply_theory_style:
        apply_theory_style(
            ax,
            xlabel=r"$|u|^{2}$",
            ylabel=r"specific energy $e$",
            title=r"Ideal-gas specific internal energy (total)",
            subtitle=r"ref $p=\rho=1$, $\gamma=1.4$ ($e_{\mathrm{int}}=2.5$)",
        )
    else:
        ax.set_xlabel(r"|u|^2")
        ax.set_ylabel(r"e")
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
