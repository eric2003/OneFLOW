#!/usr/bin/env python3
"""Roe dissipation wave coefficients xi1, xi2 vs |vn| (fixed a2,a3,c)."""
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
            xs.append(float(row["abs_vn"]))
            y1.append(float(row["xi1"]))
            y2.append(float(row["xi2"]))
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    ax.plot(xs, y1, lw=2.0, label=r"$\xi_1$")
    ax.plot(xs, y2, lw=2.0, label=r"$\xi_2$ (independent of $|v_n|$)")
    ax.axhline(0.0, color="0.6", ls=":", lw=1)
    if apply_theory_style:
        apply_theory_style(
            ax,
            xlabel=r"$a_1=|v_n|$",
            ylabel=r"wave coefficients",
            title=r"Roe dissipation wave coefficients",
            subtitle=r"ref $a_2=3$, $a_3=0.5$, $c=2$: "
                     r"$\xi_1=(2a_1-a_2-a_3)/(2c^{2})$, $\xi_2=(a_2-a_3)/(2c)$",
        )
    else:
        ax.set_xlabel(r"$a_1=|v_n|$")
        ax.set_ylabel("wave coefficients")
        ax.set_title("Roe dissipation wave coefficients")
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
