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
            xs.append(float(row["d"]))
            ys.append(float(row["omega_w"]))
    fig, ax = plt.subplots(figsize=(7.0, 4.2))
    ax.plot(xs, ys, lw=2.0, color="C4", label='$\\omega_w(d)$')
    ax.set_yscale("log")
    apply_theory_style(
        ax,
        xlabel='$d$ (wall distance)',
        ylabel='$\\omega_w$',
        title='Menter wall $\\omega_w$ vs distance',
        subtitle='$\\omega_w=60\\mu_l/(\\rho\\beta_1 d^{2} Re)$  [ref $\\mu=\\rho=Re=1$, $\\beta_1=0.075$]',
    )
    save_figure(fig, args.out)

if __name__ == "__main__":
    main()
