#!/usr/bin/env python3
"""Face Mach number vs relative normal velocity (fixed sound speed)."""
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
            xs.append(float(row["vn_rel"]))
            ys.append(float(row["mach"]))
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    ax.plot(xs, ys, lw=2.0, label=r"$M=v_n/c$")
    ax.axhline(0.0, color="0.5", ls=":", lw=1)
    if apply_theory_style:
        apply_theory_style(
            ax,
            xlabel=r"$v_n$ (relative normal velocity)",
            ylabel=r"Mach $M$",
            title=r"Face Mach number",
            subtitle=r"ref $c=1.5$: $M=v_n/c$",
        )
    else:
        ax.set_xlabel(r"vn")
        ax.set_ylabel(r"M")
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
