#!/usr/bin/env python3
"""Run one independent Python STL-to-Gmsh reference Case."""
from __future__ import annotations
import argparse, json, shutil, subprocess, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[3]
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
    results.mkdir(parents=True, exist_ok=True)
    geometry = convert(source, sidecar, geo, tolerance=1.0e-10)
    subprocess.run([gmsh, str(geo), "-2", "-format", "msh4", "-o", str(msh)],
                   check=True, cwd=case)
    mesh = inspect_msh4(msh)
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
    report = {"case": case.name, "geometry": "results/model.geo",
              "mesh": "results/model.msh", "geometry_report": geometry,
              "mesh_summary": mesh, "expected": expected, "validated": True}
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
