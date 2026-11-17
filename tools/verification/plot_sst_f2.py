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
            xs.append(float(row["arg2"]))
            ys.append(float(row["F2"]))
    fig, ax = plt.subplots(figsize=(7.0, 4.2))
    ax.plot(xs, ys, lw=2.0, color="C1", label='$F_2$')
    apply_theory_style(
        ax,
        xlabel='$\\mathrm{arg}_2$',
        ylabel='$F_2$',
        title='SST limiter $F_2$',
        subtitle='$F_2=\\tanh(\\mathrm{arg}_2^{2})$',
    )
    save_figure(fig, args.out)

if __name__ == "__main__":
    main()
