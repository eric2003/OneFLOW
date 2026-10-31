#!/usr/bin/env python3
"""Plot SA wall-destruction function fw(r) from theory CSV.

Requires matplotlib. Output is a static PNG (Agg), suitable for docs.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--csv",
        type=Path,
        default=Path(__file__).resolve().parent / "data" / "sa_fw_r_theory.csv",
    )
    ap.add_argument(
        "--out",
        type=Path,
        default=Path(__file__).resolve().parent / "data" / "sa_fw_r_theory.png",
    )
    args = ap.parse_args()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    r_vals, fw_vals = [], []
    with args.csv.open(newline="", encoding="ascii") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            r_vals.append(float(row["r"]))
            fw_vals.append(float(row["fw"]))

    fig, ax = plt.subplots(figsize=(6.0, 4.0))
    ax.plot(r_vals, fw_vals, color="C0", lw=2.0, label=r"$f_w(r)$")
    ax.axhline(1.0, color="0.5", ls="--", lw=1.0, label=r"$f_w=1$ (log layer)")
    ax.axvline(1.0, color="0.5", ls=":", lw=1.0, label=r"$r=1$")
    ax.set_xlabel(r"$r = \min(\tilde\nu/(\tilde S \kappa^2 d^2), 10)$")
    ax.set_ylabel(r"$f_w$")
    ax.set_title("Spalart–Allmaras wall destruction function")
    ax.set_xlim(0.0, 10.0)
    ax.set_ylim(0.0, None)
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best")
    fig.tight_layout()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=120)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
