#!/usr/bin/env python3
"""Roe dissipation coupling terms vs enthalpy jump dh."""
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
    xs, y1, y2 = [], [], []
    with args.csv.open(newline="", encoding="ascii") as fh:
        for row in csv.DictReader(fh):
            xs.append(float(row["dh"]))
            y1.append(float(row["term1"]))
            y2.append(float(row["term2"]))
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    ax.plot(xs, y1, lw=2.0, label=r"term$_1$")
    ax.plot(xs, y2, lw=2.0, label=r"term$_2$")
    ax.axhline(0.0, color="0.6", ls=":", lw=1)
    ax.axvline(0.0, color="0.6", ls=":", lw=1)
    if apply_theory_style:
        apply_theory_style(
            ax,
            xlabel=r"$\Delta H$ (total-enthalpy jump)",
            ylabel=r"coupling terms",
            title=r"Roe dissipation coupling terms",
            subtitle=r"ref $d_c=2$, $c^{2}=4$, $\xi_1=0.5$, $\xi_2=0.25$: "
                     r"$\mathrm{term}_1=\Delta H\xi_1+d_c\xi_2$",
        )
    else:
        ax.set_xlabel(r"$\Delta H$")
        ax.set_ylabel("coupling terms")
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
