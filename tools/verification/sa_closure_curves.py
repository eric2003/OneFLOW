#!/usr/bin/env python3
"""Generate intermediate SA closure curves from published closed forms.

This script evaluates theoretical expressions (not OneFLOW production C++).
C++ gtest compares production SaG/SaFw against the same closed forms.

Proves: algebraic agreement of pure APIs with the published fw(r) formula.
Does not prove: mesh discretization, wall distance, or full-case CFD accuracy.

Constants (Spalart-Allmaras / NASA TMR defaults):
  cw2 = 0.3, cw3 = 2.0
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def sa_g(r: float, cw2: float = 0.3) -> float:
    return r + cw2 * (r**6 - r)


def sa_fw(r: float, cw2: float = 0.3, cw3: float = 2.0) -> float:
    g = sa_g(r, cw2)
    cw36 = cw3**6
    g6 = g**6
    return g * ((1.0 + cw36) / (g6 + cw36)) ** (1.0 / 6.0)


def sa_fv1(chi: float, cv1: float = 7.1) -> float:
    chi3 = chi**3
    cv13 = cv1**3
    return chi3 / (chi3 + cv13)


def write_fw_csv(path: Path) -> None:
    rs = [0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 5.0, 10.0]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="ascii") as fh:
        w = csv.writer(fh)
        w.writerow(["r", "g", "fw"])
        for r in rs:
            w.writerow([f"{r:.12g}", f"{sa_g(r):.16e}", f"{sa_fw(r):.16e}"])


def write_fv1_csv(path: Path) -> None:
    chis = [0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 50.0, 100.0]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="ascii") as fh:
        w = csv.writer(fh)
        w.writerow(["chi", "fv1"])
        for chi in chis:
            w.writerow([f"{chi:.12g}", f"{sa_fv1(chi):.16e}"])


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--out-dir",
        type=Path,
        default=Path(__file__).resolve().parent / "data",
    )
    args = ap.parse_args()
    write_fw_csv(args.out_dir / "sa_fw_r_theory.csv")
    write_fv1_csv(args.out_dir / "sa_fv1_chi_theory.csv")
    print(f"wrote curves under {args.out_dir}")


if __name__ == "__main__":
    main()
