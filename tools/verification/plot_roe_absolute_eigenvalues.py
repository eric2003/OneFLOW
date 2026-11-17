#!/usr/bin/env python3
"""Roe absolute eigenvalues vs vn (fixed c=1)."""
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
    xs, a1, a2, a3 = [], [], [], []
    with args.csv.open(newline="", encoding="ascii") as fh:
        for row in csv.DictReader(fh):
            xs.append(float(row["vn"]))
            a1.append(float(row["abs_vn"]))
            a2.append(float(row["abs_vn_plus_c"]))
            a3.append(float(row["abs_vn_minus_c"]))
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    ax.plot(xs, a1, lw=2.0, label=r"$|v_n|$")
    ax.plot(xs, a2, lw=2.0, label=r"$|v_n+c|$")
    ax.plot(xs, a3, lw=2.0, label=r"$|v_n-c|$")
    ax.axvline(0.0, color="0.6", ls=":", lw=1)
    ax.axvline(-1.0, color="0.7", ls=":", lw=1)
    ax.axvline(1.0, color="0.7", ls=":", lw=1)
    if apply_theory_style:
        apply_theory_style(
            ax,
            xlabel=r"$v_n$ (relative face normal velocity)",
            ylabel=r"absolute eigenvalues",
            title=r"Roe absolute eigenvalues (fixed $c=1$)",
            subtitle=r"$a_1=|v_n|,\ a_2=|v_n+c|,\ a_3=|v_n-c|$",
        )
    else:
        ax.set_xlabel(r"$v_n$")
        ax.set_ylabel("absolute eigenvalues")
        ax.set_title(r"Roe absolute eigenvalues (fixed $c=1$)")
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
