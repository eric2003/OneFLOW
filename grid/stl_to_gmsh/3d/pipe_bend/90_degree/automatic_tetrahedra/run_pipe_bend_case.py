#!/usr/bin/env python3
"""Generate the 90-degree pipe-bend STL and its Gmsh tetrahedral volume mesh."""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case-directory", type=Path, required=True)
    parser.add_argument("--gmsh", required=True)
    parser.add_argument("--mesh-size", type=float, default=0.20)
    args = parser.parse_args()
    case = args.case_directory.resolve()
    input_dir = case / "input"
    results_dir = case / "results"
    input_dir.mkdir(parents=True, exist_ok=True)
    results_dir.mkdir(parents=True, exist_ok=True)
    stl = input_dir / "pipe_bend_90_degree.stl"
    msh = results_dir / "pipe_bend_90_degree.msh"
    generator = ROOT / "grid/stl_to_gmsh/3d/pipe_bend/90_degree/geometry/python/generate_pipe_bend.py"
    mesher = ROOT / "grid/stl_to_gmsh/generate_volume_mesh.py"
    commands = (
        [sys.executable, str(generator), "--output", str(stl)],
        [sys.executable, str(mesher), "--stl", str(stl), "--output", str(msh),
         "--gmsh", args.gmsh, "--mesh-size", str(args.mesh_size),
         "--boundary-name", "wall", "--volume-name", "fluid_domain"],
    )
    for command in commands:
        print("Running:", subprocess.list2cmdline(command), flush=True)
        subprocess.run(command, cwd=ROOT, check=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
