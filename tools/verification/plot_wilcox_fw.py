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
    with args.csv.open(newline="", encoding="ascii") as fh:
        reader = csv.DictReader(fh)
        field = [c for c in reader.fieldnames if c != "xw"][0]
        xs, ys = [], []
        for row in reader:
            xs.append(float(row["xw"]))
            ys.append(float(row[field]))
    fig, ax = plt.subplots(figsize=(7.0, 4.2))
    if "1998" in field:
        label = r"$f_w^{1998}$"
        title = r"Wilcox $f_w^{1998}(x_w)$"
        sub = r"$f_w=(1+70 x_w)/(1+80 x_w),\ \lim=7/8$"
    else:
        label = r"$f_w^{2006}$"
        title = r"Wilcox $f_w^{2006}(x_w)$"
        sub = r"$f_w=(1+85 x_w)/(1+100 x_w),\ \lim=0.85$"
    ax.plot(xs, ys, lw=2.0, label=label)
    apply_theory_style(ax, xlabel=r"$x_w$", ylabel=label, title=title, subtitle=sub)
    save_figure(fig, args.out)

if __name__ == "__main__":
    main()
