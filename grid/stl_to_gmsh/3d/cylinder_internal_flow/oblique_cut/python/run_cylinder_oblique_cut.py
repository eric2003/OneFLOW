#!/usr/bin/env python3
"""Run the registered blind-input oblique-inlet cylinder cutting case."""
from __future__ import annotations
import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[6]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from grid.stl_import.cut_surface_cli import main as cut_surface_main
from grid.stl_to_gmsh.generate_volume_mesh import discover_gmsh, generate

DEFAULT_CASE = Path(__file__).resolve().parent
DEFAULT_OUTPUT = Path.home() / ".oneflow/mesh-generation/cylinder_oblique_cut"


def run_case(case_directory: Path, output_directory: Path, gmsh: str) -> dict:
    case, output = case_directory.expanduser().resolve(), output_directory.expanduser().resolve()
    try:
        output.relative_to(ROOT)
    except ValueError:
        pass
    else:
        raise ValueError(f"Generated files must stay outside the source repository: {output}")
    # Prefer assets staged inside the selected implementation workspace.
    # The parent-level location remains a source-tree compatibility fallback.
    config_candidates = (
        case / "script/cutting.json",
        case.parent / "script/cutting.json",
        DEFAULT_CASE.parent / "script/cutting.json",
    )
    config_path = next((path for path in config_candidates if path.is_file()), config_candidates[0])
    source = case / "input/model.stl"
    # The native and Python implementations share one blind STL fixture at
    # the workflow level. A GUI workspace may stage only the selected
    # implementation directory, so resolve the shared source from the
    # repository case when it is not present in the staged Python workspace.
    if not source.is_file():
        repository_source = (
            DEFAULT_CASE.parent / "cpp" / "input" / "model.stl"
        ).resolve()
        if repository_source.is_file():
            source = repository_source
    missing = [
        str(path) for path in (config_path, source) if not path.is_file()
    ]
    if missing:
        raise FileNotFoundError(
            "Incomplete oblique-cut case; missing required assets: "
            + ", ".join(missing)
        )
    config = json.loads(config_path.read_text(encoding="utf-8"))
    planes = config["cutting_planes"]
    inlet, outlet = planes["inlet"], planes["outlet"]
    cut_surface_main([
        "--input", str(source), "--output-directory", str(output),
        "--inlet-point", *(str(v) for v in inlet["point"]),
        "--inlet-normal", *(str(v) for v in inlet["normal"]),
        "--outlet-point", *(str(v) for v in outlet["point"]),
        "--outlet-normal", *(str(v) for v in outlet["normal"]),
        "--inlet-keep", "positive" if inlet.get("keep_positive", True) else "negative",
        "--outlet-keep", "positive" if outlet.get("keep_positive", False) else "negative",
        "--tolerance", str(config.get("tolerance", 1.0e-8)),
    ])
    mesh_config = config["volume_mesh"]
    report = generate(
        output / "closed_surface.stl",
        output / "results/cylinder_oblique_cut.msh",
        gmsh,
        mesh_size_min=float(mesh_config["mesh_size_min"]),
        mesh_size_max=float(mesh_config["mesh_size_max"]),
    )
    expected = {(g["dimension"], g["name"]) for g in mesh_config["expected_physical_groups"]}
    actual = {(g["dimension"], g["name"]) for g in report["mesh_summary"]["physical_groups"]}
    if actual != expected:
        raise RuntimeError(f"Physical Groups {sorted(actual)} do not match {sorted(expected)}")
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case-directory", type=Path, default=DEFAULT_CASE)
    parser.add_argument("--output-directory", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--gmsh", default=discover_gmsh())
    args = parser.parse_args()
    if not args.gmsh:
        parser.error("Gmsh was not found; install it or pass --gmsh")
    try:
        report = run_case(args.case_directory, args.output_directory, args.gmsh)
    except (OSError, ValueError, KeyError, RuntimeError) as error:
        parser.exit(2, f"Oblique cylinder cutting failed: {error}\n")
    print(f"PASS: oblique cylinder cut and tetrahedral volume mesh generated: {report['mesh']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
