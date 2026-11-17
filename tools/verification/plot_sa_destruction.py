#!/usr/bin/env python3
"""SA destruction term vs nu_tilde/d^2 measure."""
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
            xs.append(float(row["nu_tilde_rs"]))
            ys.append(float(row["destruction"]))
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    ax.plot(xs, ys, lw=2.0, label=r"$D=(c_{w1}f_w-c_{b1}f_{t2})\tilde{\nu}_{rs}$")
    if apply_theory_style:
        apply_theory_style(
            ax,
            xlabel=r"$\tilde{\nu}_{rs}$",
            ylabel=r"destruction $D$",
            title=r"SA destruction term",
            subtitle=r"ref $c_{w1}k=3$, $f_w=0.5$, $c_{b1}=0.1355$, $f_{t2}=0.2$",
        )
    else:
        ax.set_xlabel(r"nu_rs")
        ax.set_ylabel(r"D")
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
