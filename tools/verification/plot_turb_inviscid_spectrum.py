#!/usr/bin/env python3
from __future__ import annotations
import argparse, csv, sys
from pathlib import Path

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    tools = Path(__file__).resolve().parent
    if str(tools) not in sys.path:
        sys.path.insert(0, str(tools))
    from plot_style import apply_theory_style, save_figure
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    xs, ys = [], []
    with args.csv.open(newline="", encoding="ascii") as fh:
        for row in csv.DictReader(fh):
            xs.append(float(row["abs_vn"]))
            ys.append(float(row["lambda_inv"]))
    fig, ax = plt.subplots(figsize=(7.0, 4.2))
    ax.plot(xs, ys, lw=2.0, color="C0", label='$\\lambda_{\\mathrm{inv}}$')
    apply_theory_style(
        ax,
        xlabel='$|v_n|$',
        ylabel='$\\lambda_{\\mathrm{inv}}$',
        title='Turbulence inviscid face spectral radius',
        subtitle='$\\lambda_{\\mathrm{inv}}=\\frac{A}{2}|v_n|$  [ref $A=2$]',
    )
    save_figure(fig, args.out)

if __name__ == "__main__":
    main()
