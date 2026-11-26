#!/usr/bin/env python3
"""Generate a parametric 90-degree pipe-bend STL and its tetrahedral MSH4 mesh."""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

REPOSITORY_ROOT = Path(__file__).resolve().parents[6]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case-directory", type=Path, required=True)
    parser.add_argument("--gmsh", required=True)
    parser.add_argument("--mesh-size", type=float, default=0.20)
    args = parser.parse_args()

    case_directory = args.case_directory.expanduser().resolve()
    input_directory = case_directory / "input"
    results_directory = case_directory / "results"
    input_directory.mkdir(parents=True, exist_ok=True)
    results_directory.mkdir(parents=True, exist_ok=True)
    stl = input_directory / "pipe_bend_90_degree.stl"
    msh = results_directory / "pipe_bend_90_degree.msh"
    geometry_generator = (
        REPOSITORY_ROOT
        / "grid/stl_to_gmsh/3d/pipe_bend/90_degree/geometry/python/generate_pipe_bend.py"
    )
    volume_mesher = REPOSITORY_ROOT / "grid/stl_to_gmsh/generate_volume_mesh.py"

    commands = (
        [sys.executable, str(geometry_generator), "--output", str(stl)],
        [
            sys.executable, str(volume_mesher),
            "--stl", str(stl), "--output", str(msh),
            "--gmsh", args.gmsh, "--mesh-size", str(args.mesh_size),
            "--boundary-name", "wall", "--volume-name", "fluid_domain",
        ],
    )
    for command in commands:
        print("Running:", subprocess.list2cmdline(command), flush=True)
        subprocess.run(command, cwd=REPOSITORY_ROOT, check=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
