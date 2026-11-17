#!/usr/bin/env python3
"""Harten entropy-fix floors vs |vn| (theory reference state)."""
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
    xs, yv, yp, ym = [], [], [], []
    with args.csv.open(newline="", encoding="ascii") as fh:
        for row in csv.DictReader(fh):
            xs.append(float(row["abs_vn"]))
            yv.append(float(row["fixed_vn"]))
            yp.append(float(row["fixed_vn_plus_c"]))
            ym.append(float(row["fixed_vn_minus_c"]))
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    ax.plot(xs, yv, lw=2.0, label=r"fixed $|v_n|$")
    ax.plot(xs, yp, lw=2.0, label=r"fixed $|v_n+c|$")
    ax.plot(xs, ym, lw=2.0, label=r"fixed $|v_n-c|$")
    ax.axhline(1.0, color="0.5", ls=":", lw=1)
    ax.axhline(0.5, color="0.5", ls=":", lw=1)
    if apply_theory_style:
        apply_theory_style(
            ax,
            xlabel=r"$|v_n|$ (input)",
            ylabel=r"fixed $|\lambda|$",
            title=r"Harten entropy fix (eigenvalue floor)",
            subtitle=r"ref: $|v_n+c|=2$, $|v_n-c|=1.5$, $c_1=0.5$, $c_2=0.25$; $m=\lambda_{\max} c_i$",
        )
    else:
        ax.set_xlabel(r"$|v_n|$ (input)")
        ax.set_ylabel(r"fixed $|\lambda|$")
        ax.set_title("Harten entropy fix (eigenvalue floor)\n"
                     r"ref: $|v_n+c|=2$, $|v_n-c|=1.5$, $c_1=0.5$, $c_2=0.25$")
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
