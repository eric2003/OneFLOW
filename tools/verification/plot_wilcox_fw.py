#!/usr/bin/env python3
"""Plot Wilcox fw1998 / fw2006 theory curves."""
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
        label, title, sub = "fw1998", "Wilcox fw1998(xw)", "fw=(1+70*xw)/(1+80*xw); limit 7/8"
    else:
        label, title, sub = "fw2006", "Wilcox fw2006(xw)", "fw=(1+85*xw)/(1+100*xw); limit 0.85"
    ax.plot(xs, ys, lw=2.0, label=label)
    apply_theory_style(
        ax,
        xlabel="xw (vortex-stretching parameter)",
        ylabel=label,
        title=title,
        subtitle=sub,
        legend=True,
    )
    save_figure(fig, args.out)
    print(f"wrote {args.out}")

if __name__ == "__main__":
    main()
