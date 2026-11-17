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
    s, l1, l2, l3 = [], [], [], []
    with args.csv.open(newline="", encoding="ascii") as fh:
        for row in csv.DictReader(fh):
            s.append(float(row["sign"]))
            l1.append(float(row["lmd1"]))
            l2.append(float(row["lmd2"]))
            l3.append(float(row["lmd3"]))
    fig, ax = plt.subplots(figsize=(7.0, 4.2))
    ax.plot(s, l1, "o-", lw=2, label=r"$\lambda_1$ ($v_n$)")
    ax.plot(s, l2, "s-", lw=2, label=r"$\lambda_2$ ($v_n+c$)")
    ax.plot(s, l3, "^-", lw=2, label=r"$\lambda_3$ ($v_n-c$)")
    apply_theory_style(
        ax,
        xlabel=r"matrix sign $s$",
        ylabel=r"split $\lambda$",
        title=r"LUSGS face eigenvalue split",
        subtitle=r"$\lambda=\frac{1}{2}(\lambda^{\mathrm{phys}}+s|\lambda|_f)$  [ref $v_n=1$, $c=2$, $|\lambda|_f=0.4$]",
    )
    save_figure(fig, args.out)

if __name__ == "__main__":
    main()
