#!/usr/bin/env python3
"""Generate a 3-D tetrahedral volume mesh directly from a closed STL shell."""
from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from grid.stl_to_gmsh.converter import Stl2DError, convert, inspect_msh4  # noqa: E402


def generate(
    stl: Path,
    output: Path,
    gmsh: str,
    mesh_size_min: float,
    mesh_size_max: float,
    boundary_name: str = "wall",
    volume_name: str = "fluid_domain",
    surface_angle_degrees: float = 40.0,
    curve_angle_degrees: float = 180.0,
) -> dict:
    stl = stl.expanduser().resolve()
    output = output.expanduser().resolve()
    if not stl.is_file():
        raise FileNotFoundError(f"STL input does not exist: {stl}")
    if output.suffix.lower() != ".msh":
        raise ValueError("output path must use the .msh extension")
    if output == stl:
        raise ValueError("output mesh must not overwrite the input STL")
    if mesh_size_min <= 0 or mesh_size_max <= 0:
        raise ValueError("mesh sizes must be positive")
    if mesh_size_min > mesh_size_max:
        raise ValueError("--mesh-size-min must not exceed --mesh-size-max")

    output.parent.mkdir(parents=True, exist_ok=True)
    geo = output.with_suffix(".geo")
    report_path = output.with_suffix(".report.json")
    config = {
        "volume_mesh": "closed_surface_automatic_tetrahedra",
        "mesh_dimension": 3,
        "mesh_size_min": mesh_size_min,
        "mesh_size_max": mesh_size_max,
        "surface_name": volume_name,
        "boundary_name": boundary_name,
        "surface_angle_degrees": surface_angle_degrees,
        "curve_angle_degrees": curve_angle_degrees,
    }
    with tempfile.TemporaryDirectory(prefix="oneflow-stl-to-gmsh-") as temporary:
        sidecar = Path(temporary) / "layout.json"
        sidecar.write_text(json.dumps(config), encoding="utf-8")
        geometry_report = convert(stl, sidecar, geo, tolerance=1.0e-10)

    subprocess.run(
        [gmsh, str(geo), "-3", "-format", "msh4", "-o", str(output)],
        check=True,
        cwd=output.parent,
    )
    mesh = inspect_msh4(output)
    if mesh["tetrahedron_elements"] <= 0:
        raise RuntimeError("Gmsh output contains no TETRA_4 volume elements")
    if mesh["triangle_elements"] <= 0:
        raise RuntimeError("Gmsh output contains no TRI_3 boundary elements")
    groups = {(item["dimension"], item["name"]) for item in mesh["physical_groups"]}
    expected_groups = {(2, boundary_name), (3, volume_name)}
    if groups != expected_groups:
        raise RuntimeError(
            f"Physical Groups {sorted(groups)} do not match {sorted(expected_groups)}"
        )

    report = {
        "schema": "oneflow.stl-to-gmsh-volume-report",
        "version": 1,
        "input_stl": str(stl),
        "geometry": str(geo),
        "mesh": str(output),
        "mesh_size_min": mesh_size_min,
        "mesh_size_max": mesh_size_max,
        "geometry_validation": geometry_report,
        "mesh_summary": mesh,
        "validated": True,
    }
    report_path.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    print(json.dumps(report, indent=2))
    return report


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Generate a tetrahedral volume mesh from a closed, genus-zero STL shell."
    )
    parser.add_argument("--stl", type=Path, required=True, help="input closed STL surface")
    parser.add_argument("--output", type=Path, required=True, help="output .msh file")
    parser.add_argument("--mesh-size-min", type=float, required=True)
    parser.add_argument("--mesh-size-max", type=float, required=True)
    parser.add_argument("--gmsh", default=shutil.which("gmsh"), help="Gmsh executable")
    parser.add_argument("--boundary-name", default="wall")
    parser.add_argument("--volume-name", default="fluid_domain")
    parser.add_argument("--surface-angle", type=float, default=40.0)
    parser.add_argument("--curve-angle", type=float, default=180.0)
    args = parser.parse_args()
    if not args.gmsh:
        parser.error("Gmsh was not found; install it or pass --gmsh")
    try:
        generate(
            args.stl, args.output, args.gmsh,
            args.mesh_size_min, args.mesh_size_max,
            args.boundary_name, args.volume_name,
            args.surface_angle, args.curve_angle,
        )
    except (OSError, ValueError, Stl2DError, subprocess.CalledProcessError, RuntimeError) as error:
        parser.exit(2, f"STL volume mesh generation failed: {error}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
