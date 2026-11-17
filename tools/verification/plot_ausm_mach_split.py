#!/usr/bin/env python3
"""AUSM FMSplit4 / FPSplit5 vs Mach (theory, ipn=+1)."""
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
    xs, y4, y5 = [], [], []
    with args.csv.open(newline="", encoding="ascii") as fh:
        for row in csv.DictReader(fh):
            xs.append(float(row["mach"]))
            y4.append(float(row["fmsplit4_plus"]))
            y5.append(float(row["fpsplit5_plus"]))
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    ax.plot(xs, y4, lw=2.0, label=r"$M_4=\mathrm{FMSplit4}(M,\beta,+1)$")
    ax.plot(xs, y5, lw=2.0, label=r"$P_5=\mathrm{FPSplit5}(M,\alpha,+1)$")
    ax.axvline(-1.0, color="0.6", ls=":", lw=1)
    ax.axvline(1.0, color="0.6", ls=":", lw=1)
    ax.axhline(0.0, color="0.6", ls="--", lw=1)
    if apply_theory_style:
        apply_theory_style(
            ax,
            xlabel=r"Mach number $M$",
            ylabel=r"split function value",
            title=r"AUSM Mach / pressure split polynomials",
            subtitle=r"$M_4=M_1$ ($|M|>=1$) else $M_2(1-16\beta M_2)$; $\beta=1/8$, $\alpha=3/16$",
        )
    else:
        ax.set_xlabel(r"Mach number $M$")
        ax.set_ylabel("split function value")
        ax.set_title("AUSM Mach / pressure split polynomials\n"
                     r"$\beta=1/8$, $\alpha=3/16$; sonic lines at $|M|=1$")
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
