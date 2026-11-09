#!/usr/bin/env python3
"""Run one independent Python STL-to-Gmsh reference Case."""
from __future__ import annotations
import argparse, json, shutil, subprocess, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from grid.stl.stl_to_gmsh import convert, inspect_msh4  # noqa: E402

def run_case(case_directory: Path, gmsh: str) -> dict:
    case = case_directory.expanduser().resolve()
    source, sidecar = case / "input/model.stl", case / "script/layout.json"
    expected_path = case / "expected.json"
    results = case / "results"
    geo, msh = results / "model.geo", results / "model.msh"
    for path in (source, sidecar, expected_path):
        if not path.is_file():
            raise FileNotFoundError(f"Required Case input does not exist: {path}")
    expected = json.loads(expected_path.read_text(encoding="utf-8"))
    layout = json.loads(sidecar.read_text(encoding="utf-8"))
    results.mkdir(parents=True, exist_ok=True)
    geometry = convert(source, sidecar, geo, tolerance=1.0e-10)
    subprocess.run([gmsh, str(geo), "-2", "-format", "msh4", "-o", str(msh)],
                   check=True, cwd=case)
    mesh = inspect_msh4(msh)
    domain = layout.get("domain")
    if domain is not None:
        bounds = mesh["coordinate_bounds"]
        tolerance = 1.0e-10
        for axis, lower, upper in (
            (0, "xmin", "xmax"), (1, "ymin", "ymax")
        ):
            if (
                abs(bounds["minimum"][axis] - float(domain[lower])) > tolerance
                or abs(bounds["maximum"][axis] - float(domain[upper])) > tolerance
            ):
                raise RuntimeError(
                    f"Generated mesh bounds on axis {axis} do not match the computational domain"
                )
        if abs(bounds["maximum"][2] - bounds["minimum"][2]) > tolerance:
            raise RuntimeError("Generated external-flow mesh is not planar in XY")
        quality_limits = layout.get("quality", {})
        measured_quality = mesh["triangle_quality"]
        minimum_quality = float(quality_limits.get("minimum_triangle_quality", 0.01))
        minimum_angle = float(quality_limits.get("minimum_angle_degrees", 1.0))
        if measured_quality["minimum_quality"] < minimum_quality:
            raise RuntimeError(
                f"Minimum triangle quality {measured_quality['minimum_quality']:.6g} "
                f"is below {minimum_quality:.6g}"
            )
        if measured_quality["minimum_angle_degrees"] < minimum_angle:
            raise RuntimeError(
                f"Minimum triangle angle {measured_quality['minimum_angle_degrees']:.6g} "
                f"is below {minimum_angle:.6g} degrees"
            )
    observed = {key: mesh[key] for key in (
        "nodes", "line_elements", "triangle_elements", "quadrilateral_elements")}
    groups = {(g["dimension"], g["name"]) for g in mesh["physical_groups"]}
    wanted_groups = {(g["dimension"], g["name"]) for g in expected["physical_groups"]}
    expected_mesh = expected.get("mesh", {})
    for key, value in expected_mesh.items():
        if observed.get(key) != value:
            raise RuntimeError(
                f"Mesh count {key}={observed.get(key)} does not match {value}"
            )
    for key, minimum in expected.get("mesh_minimum", {}).items():
        if observed.get(key, 0) < minimum:
            raise RuntimeError(
                f"Mesh count {key}={observed.get(key, 0)} is below minimum {minimum}"
            )
    if groups != wanted_groups:
        raise RuntimeError(f"Physical Groups {sorted(groups)} do not match {sorted(wanted_groups)}")
    report = {"schema": "oneflow.stl-to-gmsh-validation-report", "version": 1,
              "case": case.name, "geometry": "results/model.geo",
              "mesh": "results/model.msh",
              "contract": {"domain": layout.get("domain"),
                           "surface_name": layout.get("surface_name"),
                           "body_boundary_name": layout.get("body_boundary_name"),
                           "quality_limits": quality_limits},
              "geometry_report": geometry, "mesh_summary": mesh,
              "expected": expected, "validated": True}
    (results / "report.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n",
        encoding="utf-8", newline="\n")
    print(json.dumps(report, indent=2))
    return report

def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case-directory", type=Path, required=True)
    parser.add_argument("--gmsh", default=shutil.which("gmsh"))
    args = parser.parse_args()
    if not args.gmsh:
        parser.error("Gmsh executable was not found; install Gmsh or pass --gmsh")
    try:
        run_case(args.case_directory, args.gmsh)
    except (OSError, ValueError, KeyError, subprocess.CalledProcessError, RuntimeError) as error:
        parser.exit(2, f"STL-to-Gmsh Python Case failed: {error}\n")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
