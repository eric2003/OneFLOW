#!/usr/bin/env python3
"""Reject MPI measurements where one or more ranks do no solver work."""

import argparse
import csv
import glob
import math
import pathlib
import re
import sys
from typing import Dict


ACTIVE_CATEGORIES = {
    "gradient",
    "limiter",
    "reconstruction",
    "boundary_reconstruction",
    "host_pack",
    "flux_backend",
    "geometry_connectivity_H2D",
    "H2D",
    "HIP kernel",
    "D2H",
    "viscous_turbulence",
    "rk_update",
    "residual_state_update",
    "load_q",
    "calc_time_step",
    "load_residuals",
    "update_residuals",
    "calc_lhs",
    "update_flowfield",
    "calc_boundary",
    "boundary_gamma_inner",
    "boundary_viscosity_inner",
    "boundary_bc",
    "boundary_gamma_ghost",
    "boundary_viscosity_ghost",
    "boundary_set_id",
    "boundary_prepare_data",
    "boundary_calc_face_bc",
    "boundary_update_bc",
}

RANK_SUFFIX = re.compile(r"-(\d+)\.tsv$")
MIN_ACTIVE_FRACTION = 1.0e-4


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--expected-ranks", type=int, required=True)
    parser.add_argument("--label", default="measurement")
    parser.add_argument("--pattern", required=True)
    return parser.parse_args()


def read_solver_total(path: pathlib.Path) -> float:
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        if reader.fieldnames != ["category", "elapsed_ms"]:
            raise ValueError(f"{path}: unexpected header {reader.fieldnames!r}")

        total = 0.0
        for row in reader:
            category = row["category"]
            try:
                elapsed = float(row["elapsed_ms"])
            except ValueError as error:
                raise ValueError(
                    f"{path}: invalid elapsed_ms for {category!r}"
                ) from error
            if not math.isfinite(elapsed) or elapsed < 0.0:
                raise ValueError(
                    f"{path}: non-finite or negative elapsed_ms for {category!r}"
                )
            if category in ACTIVE_CATEGORIES:
                total += elapsed
    return total


def main() -> int:
    args = parse_args()
    if args.expected_ranks <= 0:
        print("expected ranks must be positive", file=sys.stderr)
        return 2

    paths = [pathlib.Path(value) for value in sorted(glob.glob(args.pattern))]
    # The opt-in per-task diagnostic is written next to the stage breakdown as
    # "<breakdown>.tasks.tsv".  It shares the prefix but is not part of the
    # per-rank stage contract.
    paths = [path for path in paths if not path.name.endswith(".tasks.tsv")]
    by_rank = {}  # type: Dict[int, pathlib.Path]
    for path in paths:
        match = RANK_SUFFIX.search(path.name)
        if match is None:
            print(f"{path}: cannot determine MPI rank from filename", file=sys.stderr)
            return 2
        rank = int(match.group(1))
        if rank in by_rank:
            print(f"duplicate breakdown for rank {rank}", file=sys.stderr)
            return 2
        by_rank[rank] = path

    expected = set(range(args.expected_ranks))
    actual = set(by_rank)
    if actual != expected:
        missing = sorted(expected - actual)
        extra = sorted(actual - expected)
        print(
            f"RANK_UTILIZATION_FAIL label={args.label} "
            f"expected_ranks={args.expected_ranks} files={len(paths)} "
            f"missing={missing} extra={extra}",
            file=sys.stderr,
        )
        return 1

    totals = {}  # type: Dict[int, float]
    try:
        for rank, path in by_rank.items():
            totals[rank] = read_solver_total(path)
    except (OSError, ValueError) as error:
        print(
            f"RANK_UTILIZATION_FAIL label={args.label} error={error}",
            file=sys.stderr,
        )
        return 1

    maximum = max(totals.values())
    minimum_required = maximum * MIN_ACTIVE_FRACTION
    idle = [
        rank
        for rank, total in sorted(totals.items())
        if maximum <= 0.0 or total < minimum_required
    ]
    if idle:
        print(
            f"RANK_UTILIZATION_FAIL label={args.label} "
            f"expected_ranks={args.expected_ranks} "
            f"active_ranks={args.expected_ranks - len(idle)} "
            f"idle_ranks={idle} "
            f"minimum_required_ms={minimum_required:.6f} "
            f"max_solver_ms={maximum:.6f}",
            file=sys.stderr,
        )
        return 1

    print(
        f"RANK_UTILIZATION_PASS label={args.label} "
        f"expected_ranks={args.expected_ranks} active_ranks={len(totals)} "
        f"min_solver_ms={min(totals.values()):.6f} "
        f"max_solver_ms={maximum:.6f}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
