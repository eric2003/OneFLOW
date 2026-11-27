#!/usr/bin/env python3
"""Generate the OECD/NEA-Vattenfall T-junction STL and automatic TETRA_4 mesh."""
from __future__ import annotations

import argparse
import json
import shutil
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[6]
GENERATOR = ROOT / "grid/stl_to_gmsh/3d/pipe_tee/oecd_nea_vattenfall/geometry/python/generate_pipe_tee.py"


def mesh_occ_volume(geo: Path, msh: Path, gmsh: str, geometry: dict) -> dict:
    """Mesh the original OpenCASCADE BooleanUnion directly, never reconstruct it from STL."""
    from grid.stl_to_gmsh.converter import inspect_msh4

    gmsh_executable = shutil.which(gmsh) or gmsh
    msh.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        [gmsh_executable, str(geo), "-3", "-format", "msh4", "-o", str(msh)],
        cwd=msh.parent, check=True, capture_output=True, text=True,
    )
    summary = inspect_msh4(msh, geometry_contract={
        "main_diameter": geometry["main_pipe"]["inner_diameter"],
        "branch_diameter": geometry["branch_pipe"]["inner_diameter"],
        "main_inlet_length": geometry["main_pipe"]["inlet_length"],
        "main_outlet_length": geometry["main_pipe"]["outlet_length"],
        "branch_inlet_length": geometry["branch_pipe"]["inlet_length"],
    })
    expected_groups = {
        (2, "wall"), (2, "main_inflow"), (2, "branch_inflow"),
        (2, "outflow"), (3, "fluid_domain"),
    }
    actual_groups = {(item["dimension"], item["name"]) for item in summary["physical_groups"]}
    if actual_groups != expected_groups:
        raise RuntimeError(f"Direct OCC T-junction Physical Groups {sorted(actual_groups)} do not match {sorted(expected_groups)}")
    if summary["tetrahedron_elements"] <= 0 or summary["triangle_elements"] <= 0:
        raise RuntimeError("Direct OCC T-junction mesh must contain tetrahedra and boundary triangles.")
    if summary["geometry_validation"]["status"] != "PASS":
        raise RuntimeError("Direct OCC T-junction volume geometry validation did not pass.")
    return summary


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
    benchmark_path = GENERATOR.parents[2] / "automatic_tetrahedra" / "benchmark.json"
    geometry = json.loads(benchmark_path.read_text(encoding="utf-8"))["geometry"]
    mesh_occ_volume(stl.with_suffix(".geo"), msh, gmsh, geometry)
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
