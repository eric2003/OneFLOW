#!/usr/bin/env python3
"""Run and validate a configured OneFLOW mesh-generation case."""

import argparse
import json
import math
import subprocess
import sys
from pathlib import Path


def read_json(path):
    with path.open("r", encoding="utf-8") as stream:
        return json.load(stream)


def check_equal(name, actual, expected):
    if actual != expected:
        raise ValueError(f"Case check failed for {name}: expected {expected}, got {actual}")


def check_close(name, actual, expected, tolerance=1.0e-12):
    if not math.isfinite(actual) or abs(actual - expected) > tolerance:
        raise ValueError(f"Case check failed for {name}: expected {expected}, got {actual}")


def validate_summary(summary, expected, prefix):
    check_equal(f"{prefix} dimension", summary["dimension"], expected["dimension"])
    counts = summary["counts"]
    zone_count = len(summary["zones"]) if "zones" in summary else counts["zones"]
    check_equal(f"{prefix} zones", zone_count, expected["zones"])
    for key in ("nodes", "cells", "faces", "boundaryFaces"):
        check_equal(f"{prefix} {key}", counts[key], expected[key])

    internal_faces = counts["faces"] - counts["boundaryFaces"]
    check_equal(f"{prefix} internalFaces", internal_faces, expected["internalFaces"])

    bounds = summary["boundingBox"]
    for axis in range(3):
        check_close(f"{prefix} boundingBox.min[{axis}]",
                    bounds["min"][axis], expected["boundingBox"]["min"][axis])
        check_close(f"{prefix} boundingBox.max[{axis}]",
                    bounds["max"][axis], expected["boundingBox"]["max"][axis])


def run_case(case_file, workdir, executable, repository_root):
    case_file = case_file.resolve()
    workdir = workdir.resolve()
    executable = executable.resolve()
    repository_root = repository_root.resolve()

    case = read_json(case_file)
    parameters = case["parameters"]
    runtime = parameters["runtime"]
    inspection = parameters["inspection"]
    expected = parameters["expected"]

    for relative_path, label in (
        ("script/control.txt", "case control file"),
        (runtime["layout"], "grid layout file"),
    ):
        path = workdir / relative_path
        if not path.is_file():
            raise FileNotFoundError(f"{label} not found: {path}")

    exporter = repository_root / "tools" / "mesh_export" / "export_structured_grid.py"
    if not exporter.is_file():
        raise FileNotFoundError(f"Mesh exporter not found: {exporter}")
    if not executable.is_file():
        raise FileNotFoundError(f"OneFLOW executable not found: {executable}")

    result = subprocess.run([str(executable), "r", str(workdir)], cwd=executable.parent, check=False)
    if result.returncode != 0:
        raise RuntimeError(f"OneFLOW mesh case failed with exit code {result.returncode}")

    generated_grid = workdir / inspection["sourceGrid"]
    oneflow_grid = workdir / inspection["oneFlowGrid"]
    grid_info = workdir / inspection["gridInfo"]
    package_dir = workdir / inspection["packageDirectory"]
    for output_path in (generated_grid, oneflow_grid, grid_info):
        if not output_path.is_file():
            raise FileNotFoundError(f"Expected OneFLOW output was not found: {output_path}")

    validate_summary(read_json(grid_info), expected, "runtime")

    export_result = subprocess.run(
        [sys.executable, str(exporter), str(generated_grid),
         "--output", str(package_dir), "--real-bytes", str(inspection["realBytes"]),
         "--byte-order", inspection["byteOrder"]],
        check=False,
    )
    if export_result.returncode != 0:
        raise RuntimeError(f"Mesh export failed with exit code {export_result.returncode}")

    summary_path = package_dir / "summary.json"
    if not summary_path.is_file():
        raise FileNotFoundError(f"Mesh package summary not found: {summary_path}")
    validate_summary(read_json(summary_path), expected, "export")

    for relative_path in ("geometry.json", "topology.json", "boundary.json",
                          "manifest.json", "views/mesh.svg", "views/mesh.vtu"):
        output_path = package_dir / relative_path
        if not output_path.is_file():
            raise FileNotFoundError(f"Expected mesh package output was not found: {output_path}")

    print("Minimal rectangle case passed.")
    print(f"Generated grid: {generated_grid}")
    print(f"OneFLOW grid: {oneflow_grid}")
    print(f"Grid summary: {grid_info}")
    print(f"Mesh package: {package_dir}")
    print(f"SVG view: {package_dir / 'views' / 'mesh.svg'}")
    print(f"VTU view: {package_dir / 'views' / 'mesh.vtu'}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    run_parser = subparsers.add_parser("run", help="run and inspect a mesh case")
    run_parser.add_argument("--case", type=Path, required=True,
                            help="case.json from the repository")
    run_parser.add_argument("--workdir", type=Path, required=True,
                            help="writable case directory containing script/ and inputs")
    run_parser.add_argument("--oneflow", type=Path, required=True,
                            help="path to the built OneFLOW executable")
    run_parser.add_argument("--repository-root", type=Path,
                            default=Path(__file__).resolve().parents[1],
                            help="OneFLOW repository root")

    args = parser.parse_args()
    try:
        if args.command == "run":
            run_case(args.case, args.workdir, args.oneflow, args.repository_root)
    except (OSError, KeyError, ValueError, RuntimeError, json.JSONDecodeError) as error:
        print(f"Mesh workflow failed: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
