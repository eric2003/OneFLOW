#!/usr/bin/env python3
"""Generate one- and two-cell MSH 4.1 meshes from the bundled planar STL."""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path

REPOSITORY_ROOT = Path(__file__).resolve().parents[3]
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

from grid.stl.stl_to_gmsh import convert, inspect_msh4  # noqa: E402


EXAMPLES = {
    "single_quad": {
        "input": "square.stl", "x": 1, "y": 1, "nodes": 4, "lines": 4, "quads": 1
    },
    "two_quads": {
        "input": "two_square_cells.stl", "x": 2, "y": 1, "nodes": 6, "lines": 6, "quads": 2
    },
}
REQUIRED_GROUPS = {"inflow", "outflow", "symmetry", "farfield", "fluid_domain"}


def run_example(case_dir: Path, name: str, gmsh: str) -> dict:
    expected = EXAMPLES[name]
    stl = case_dir / "input" / expected["input"]
    sidecar = case_dir / "script" / f"{name}.json"
    results = case_dir / "results"
    geo = results / f"{name}.geo"
    msh = results / f"{name}.msh"
    report_path = results / f"{name}.report.json"

    conversion = convert(stl, sidecar, geo)
    conversion["geo"] = str(geo.relative_to(case_dir))
    subprocess.run(
        [gmsh, str(geo), "-2", "-format", "msh4", "-o", str(msh)],
        check=True,
        cwd=case_dir,
    )
    mesh = inspect_msh4(msh)
    names = {group["name"] for group in mesh["physical_groups"]}
    observed = {
        "nodes": mesh["nodes"],
        "lines": mesh["line_elements"],
        "triangles": mesh["triangle_elements"],
        "quads": mesh["quadrilateral_elements"],
    }
    wanted = {
        "nodes": expected["nodes"],
        "lines": expected["lines"],
        "triangles": 0,
        "quads": expected["quads"],
    }
    physical_groups = {
        (group["dimension"], group["name"]) for group in mesh["physical_groups"]
    }
    expected_groups = {(1, name) for name in REQUIRED_GROUPS - {"fluid_domain"}}
    expected_groups.add((2, "fluid_domain"))
    if observed != wanted:
        raise RuntimeError(f"{name}: generated mesh counts {observed} do not match {wanted}")
    if names != REQUIRED_GROUPS or physical_groups != expected_groups:
        raise RuntimeError(
            f"{name}: Physical Groups {sorted(physical_groups)} do not match "
            f"{sorted(expected_groups)}"
        )
    report = {
        "case": name,
        "source_stl": f"input/{expected['input']}",
        "geometry": str(geo.relative_to(case_dir)),
        "mesh": str(msh.relative_to(case_dir)),
        "conversion": conversion,
        "mesh_summary": mesh,
        "expected": wanted,
        "validated": True,
    }
    report_path.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--case-directory",
        type=Path,
        required=True,
        help="staged copy of grid/stl/stl_to_gmsh",
    )
    parser.add_argument(
        "--example",
        choices=("all", *EXAMPLES),
        default="all",
        help="generate one example or both (default: all)",
    )
    parser.add_argument(
        "--gmsh",
        default=shutil.which("gmsh"),
        help="Gmsh executable (default: search PATH)",
    )
    args = parser.parse_args()
    if not args.gmsh:
        parser.error("Gmsh executable was not found; install Gmsh or pass --gmsh")
    case_dir = args.case_directory.expanduser().resolve()
    for relative in (
        "input/square.stl", "input/two_square_cells.stl",
        "script/single_quad.json", "script/two_quads.json",
    ):
        if not (case_dir / relative).is_file():
            parser.error(f"case directory is missing {relative}: {case_dir}")
    (case_dir / "results").mkdir(parents=True, exist_ok=True)

    names = list(EXAMPLES) if args.example == "all" else [args.example]
    try:
        reports = [run_example(case_dir, name, args.gmsh) for name in names]
    except (OSError, ValueError, subprocess.CalledProcessError, RuntimeError) as exc:
        parser.exit(2, f"stl_to_gmsh case failed: {exc}\n")
    print(json.dumps(reports, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
