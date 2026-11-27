#!/usr/bin/env python3
"""Generate the OECD/NEA-Vattenfall T-junction STL and automatic TETRA_4 mesh."""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[6]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
GENERATOR = ROOT / "grid/stl_to_gmsh/3d/pipe_tee/oecd_nea_vattenfall/geometry/python/generate_pipe_tee.py"


OPERATING_MODES = {
    "benchmark": "Original OECD/NEA-Vattenfall benchmark assignment: one main inlet and one branch inlet.",
    "main_inlets": "Both ends of the straight main pipe are inlets; the vertical branch is the outlet.",
    "branch_inlet": "The vertical branch is the inlet; both ends of the straight main pipe are outlets.",
}


def configure_operating_mode(geo: Path, operating_mode: str) -> None:
    """Assign symmetric main-pipe ports and the opposite branch boundary type."""
    if operating_mode not in OPERATING_MODES:
        raise ValueError(f"Unsupported T-junction operating mode: {operating_mode}")
    if operating_mode == "benchmark":
        return
    source = geo.read_text(encoding="utf-8")
    old_groups = (
        'Physical Surface("main_inflow", 2) = {inlet[]};',
        'Physical Surface("branch_inflow", 3) = {branch_inlet[]};',
        'Physical Surface("outflow", 4) = {outlet[]};',
    )
    if not all(item in source for item in old_groups):
        raise ValueError("The generated T-junction geometry does not contain the expected three-port Physical Groups.")
    for item in old_groups:
        source = source.replace(item, "")
    if operating_mode == "main_inlets":
        groups = """Physical Surface("main_left_inflow", 2) = {inlet[]};
Physical Surface("main_right_inflow", 4) = {outlet[]};
Physical Surface("branch_outflow", 5) = {branch_inlet[]};"""
    else:
        groups = """Physical Surface("main_left_outflow", 2) = {inlet[]};
Physical Surface("main_right_outflow", 4) = {outlet[]};
Physical Surface("branch_inflow", 3) = {branch_inlet[]};"""
    anchor = 'Physical Surface("wall", 1) = {wall[]};'
    if anchor not in source:
        raise ValueError("The generated T-junction geometry is missing its wall Physical Group.")
    source = source.replace(anchor, anchor + "\n" + groups)
    geo.write_text(source, encoding="ascii")


def mesh_occ_volume(geo: Path, msh: Path, gmsh: str, geometry: dict,
                    operating_mode: str = "benchmark") -> dict:
    """Mesh the original OpenCASCADE BooleanUnion directly, never reconstruct it from STL."""
    configure_operating_mode(geo, operating_mode)
    from grid.stl_to_gmsh.converter import inspect_msh4

    gmsh_executable = shutil.which(gmsh) or gmsh
    msh.parent.mkdir(parents=True, exist_ok=True)
    volume_geo = msh.with_suffix(".geo")
    volume_geo.write_text(
        geo.read_text(encoding="utf-8") + "\nMesh.MshFileVersion = 4.1;\n",
        encoding="utf-8",
    )
    subprocess.run(
        [gmsh_executable, str(volume_geo), "-3", "-format", "msh4", "-o", str(msh)],
        cwd=msh.parent, check=True, capture_output=True, text=True,
    )
    summary = inspect_msh4(msh, geometry_contract={
        "main_diameter": geometry["main_pipe"]["inner_diameter"],
        "branch_diameter": geometry["branch_pipe"]["inner_diameter"],
        "main_inlet_length": geometry["main_pipe"]["inlet_length"],
        "main_outlet_length": geometry["main_pipe"]["outlet_length"],
        "branch_inlet_length": geometry["branch_pipe"]["inlet_length"],
    })
    if operating_mode == "benchmark":
        expected_groups = {
            (2, "wall"), (2, "main_inflow"), (2, "branch_inflow"),
            (2, "outflow"), (3, "fluid_domain"),
        }
    elif operating_mode == "main_inlets":
        expected_groups = {
            (2, "wall"), (2, "main_left_inflow"), (2, "main_right_inflow"),
            (2, "branch_outflow"), (3, "fluid_domain")
        }
    else:
        expected_groups = {
            (2, "wall"), (2, "main_left_outflow"), (2, "main_right_outflow"),
            (2, "branch_inflow"), (3, "fluid_domain")
        }
    actual_groups = {(item["dimension"], item["name"]) for item in summary["physical_groups"]}
    if actual_groups != expected_groups:
        raise RuntimeError(f"Direct OCC T-junction Physical Groups {sorted(actual_groups)} do not match {sorted(expected_groups)}")
    group_counts = {
        (item["dimension"], item["name"]): item.get("element_count", 0)
        for item in summary["physical_groups"]
    }
    empty_groups = sorted(
        f"{dimension}:{name}" for dimension, name in expected_groups
        if group_counts.get((dimension, name), 0) <= 0
    )
    if empty_groups:
        raise RuntimeError(f"Direct OCC T-junction Physical Groups are empty: {empty_groups}")
    if summary["tetrahedron_elements"] <= 0 or summary["triangle_elements"] <= 0:
        raise RuntimeError("Direct OCC T-junction mesh must contain tetrahedra and boundary triangles.")
    if summary["geometry_validation"]["status"] != "PASS":
        raise RuntimeError("Direct OCC T-junction volume geometry validation did not pass.")
    report_path = msh.with_suffix(".geometry-validation.json")
    report = {
        "schema": "oneflow.pipe-tee-occ-volume-validation",
        "version": 1,
        "geometry": "oecd_nea_vattenfall_t_junction",
        "geometry_script": str(volume_geo),
        "mesh_file": str(msh),
        "mesh_sha256": hashlib.sha256(msh.read_bytes()).hexdigest(),
        "geometry_contract": geometry,
        "mesh_summary": summary,
    }
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    summary["validation_report"] = str(report_path)
    return summary


def run_case(case_directory: Path, *, gmsh: str, mesh_size: float = 0.03,
             operating_mode: str = "benchmark", meshing_route: str = "direct_occ") -> Path:
    if meshing_route not in {"direct_occ", "stl"}:
        raise ValueError(f"Unsupported T-junction meshing route: {meshing_route}")
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
    if meshing_route == "stl":
        from grid.stl_to_gmsh.generate_volume_mesh import generate as generate_volume_mesh

        report = generate_volume_mesh(stl, msh, gmsh, mesh_size, mesh_size)
        if report.get("validated") is not True:
            raise RuntimeError("Generic STL-derived T-junction volume meshing failed.")
    else:
        mesh_occ_volume(stl.with_suffix(".geo"), msh, gmsh, geometry, operating_mode)
    return msh


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case-directory", type=Path, required=True)
    parser.add_argument("--gmsh", required=True)
    parser.add_argument("--mesh-size", type=float, default=0.03)
    parser.add_argument("--operating-mode", choices=tuple(OPERATING_MODES), default="benchmark")
    parser.add_argument("--meshing-route", choices=("direct_occ", "stl"), default="direct_occ")
    args = parser.parse_args()
    try:
        mesh = run_case(
            args.case_directory, gmsh=args.gmsh, mesh_size=args.mesh_size,
            operating_mode=args.operating_mode, meshing_route=args.meshing_route,
        )
    except (OSError, ValueError, RuntimeError, subprocess.CalledProcessError) as error:
        parser.exit(2, f"T-junction Case failed: {error}\n")
    print(f"PASS: generated {mesh}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
