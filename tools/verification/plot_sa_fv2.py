#!/usr/bin/env python3
"""Plot SA fv2(chi) theory curve (closed form, not production C++)."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--csv",
        type=Path,
        default=Path(__file__).resolve().parent / "data" / "sa_fv2_chi_theory.csv",
    )
    ap.add_argument(
        "--out",
        type=Path,
        default=Path(__file__).resolve().parent / "data" / "sa_fv2_chi_theory.png",
    )
    args = ap.parse_args()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    chi, fv2 = [], []
    with args.csv.open(newline="", encoding="ascii") as fh:
        for row in csv.DictReader(fh):
            chi.append(float(row["chi"]))
            fv2.append(float(row["fv2"]))

    fig, ax = plt.subplots(figsize=(6.0, 4.0))
    ax.plot(chi, fv2, color="C2", lw=2.0, label=r"$f_{v2}(\chi)$")
    ax.axhline(0.0, color="0.5", ls="--", lw=1.0)
    ax.axvline(7.1, color="0.5", ls=":", lw=1.0, label=r"$\chi=c_{v1}$")
    ax.set_xlabel(r"$\chi=\tilde\nu/\nu$")
    ax.set_ylabel(r"$f_{v2}$")
    ax.set_title("Spalart–Allmaras $f_{v2}$ (strain correction)")
    ax.set_xlim(0.0, 20.0)
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best")
    fig.tight_layout()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=120)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
