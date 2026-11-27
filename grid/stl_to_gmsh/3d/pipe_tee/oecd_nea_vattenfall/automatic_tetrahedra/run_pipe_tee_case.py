#!/usr/bin/env python3
"""Generate the OECD/NEA-Vattenfall T-junction STL and automatic TETRA_4 mesh."""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[6]
GENERATOR = ROOT / "grid/stl_to_gmsh/3d/pipe_tee/oecd_nea_vattenfall/geometry/python/generate_pipe_tee.py"
VOLUME_MESHER = ROOT / "grid/stl_to_gmsh/generate_volume_mesh.py"


def run_case(case_directory: Path, *, gmsh: str, mesh_size: float = 0.03) -> Path:
    case_directory = case_directory.expanduser().resolve()
    input_dir = case_directory / "input"
    results_dir = case_directory / "results"
    input_dir.mkdir(parents=True, exist_ok=True)
    results_dir.mkdir(parents=True, exist_ok=True)
    stl = input_dir / "pipe_tee_oecd_nea_vattenfall.stl"
    msh = results_dir / "pipe_tee_oecd_nea_vattenfall.msh"
    subprocess.run(
        [sys.executable, str(GENERATOR), "--output", str(stl), "--gmsh", gmsh, "--mesh-size", str(mesh_size)],
        cwd=ROOT, check=True,
    )
    subprocess.run(
        [sys.executable, str(VOLUME_MESHER), "--stl", str(stl), "--output", str(msh),
         "--gmsh", gmsh, "--mesh-size", str(mesh_size),
         "--boundary-name", "wall", "--volume-name", "fluid_domain"],
        cwd=ROOT, check=True,
    )
    return msh


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case-directory", type=Path, required=True)
    parser.add_argument("--gmsh", required=True)
    parser.add_argument("--mesh-size", type=float, default=0.03)
    args = parser.parse_args()
    try:
        mesh = run_case(args.case_directory, gmsh=args.gmsh, mesh_size=args.mesh_size)
    except (OSError, ValueError, subprocess.CalledProcessError) as error:
        parser.exit(2, f"T-junction Case failed: {error}\n")
    print(f"PASS: generated {mesh}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
