#!/usr/bin/env python3
"""Euler convective face fluxes vs relative normal velocity."""
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
    xs, fr, fx, fe = [], [], [], []
    with args.csv.open(newline="", encoding="ascii") as fh:
        for row in csv.DictReader(fh):
            xs.append(float(row["vn_rel"]))
            fr.append(float(row["flux_rho"]))
            fx.append(float(row["flux_mx"]))
            fe.append(float(row["flux_e"]))
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    ax.plot(xs, fr, lw=2.0, label=r"$\rho v_n$")
    ax.plot(xs, fx, lw=2.0, label=r"$f_{mx}$")
    ax.plot(xs, fe, lw=2.0, label=r"$f_E$")
    if apply_theory_style:
        apply_theory_style(
            ax,
            xlabel=r"$v_n$ relative",
            ylabel=r"face flux",
            title=r"Euler convective face flux",
            subtitle=r"ref $\rho=1.2$, $u=2$, $p=3$, $h=10$, $\mathbf{n}=(1,0,0)$",
        )
    else:
        ax.set_xlabel(r"vn")
        ax.set_ylabel(r"flux")
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
