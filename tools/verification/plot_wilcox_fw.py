#!/usr/bin/env python3
"""Plot Wilcox fw1998 / fw2006 theory curves."""
from __future__ import annotations
import argparse, csv
from pathlib import Path

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--variant", choices=("1998", "2006"), default=None)
    args = ap.parse_args()
    # Infer column from CSV header
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
    label = r"$f_w^{1998}$" if "1998" in field else r"$f_w^{2006}$"
    fig, ax = plt.subplots(figsize=(6.5, 4.0))
    ax.plot(xs, ys, lw=2.0, label=label)
    ax.axhline(1.0, color="0.5", ls="--", lw=1)
    ax.set_xlabel(r"$x_w$")
    ax.set_ylabel(r"$f_w$")
    ax.set_title(f"Wilcox {field}")
    ax.set_xlim(0, 5)
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best")
    fig.tight_layout()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=140)

if __name__ == "__main__":
    main()
