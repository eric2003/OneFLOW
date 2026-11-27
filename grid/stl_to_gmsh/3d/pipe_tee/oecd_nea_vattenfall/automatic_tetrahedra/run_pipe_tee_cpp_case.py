#!/usr/bin/env python3
"""Compile and run the C++17 T-junction geometry driver and direct OCC volume meshing."""
from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[6]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from run_pipe_tee_case import mesh_occ_volume

SOURCE = ROOT / "grid/stl_to_gmsh/3d/pipe_tee/oecd_nea_vattenfall/geometry/cpp/generate_pipe_tee.cpp"
EXTRACTOR = ROOT / "grid/stl_to_gmsh/3d/pipe_tee/oecd_nea_vattenfall/geometry/python/generate_pipe_tee.py"


def discover_compiler() -> str:
    requested = os.environ.get("CXX", "").strip()
    if requested:
        compiler = shutil.which(requested)
        if not compiler:
            raise RuntimeError(f"Configured CXX compiler was not found: {requested}")
        return compiler
    for name in ("c++", "g++", "clang++", "cl"):
        compiler = shutil.which(name)
        if compiler:
            return compiler
    raise RuntimeError("A C++17 compiler is required for the native geometry driver.")


def build(build_dir: Path, compiler: str | None = None) -> Path:
    compiler = compiler or discover_compiler()
    build_dir.mkdir(parents=True, exist_ok=True)
    executable = build_dir / ("pipe_tee_generator.exe" if os.name == "nt" else "pipe_tee_generator")
    if Path(compiler).name.lower() in {"cl", "cl.exe"}:
        command = [compiler, "/nologo", "/std:c++17", "/EHsc", "/O2", str(SOURCE), f"/Fe:{executable}"]
    else:
        command = [compiler, "-std=c++17", "-O2", str(SOURCE), "-o", str(executable)]
    subprocess.run(command, cwd=build_dir, check=True)
    return executable


def run_case(case_directory: Path, *, gmsh: str, mesh_size: float = 0.03,
             compiler: str | None = None, operating_mode: str = "benchmark",
             meshing_route: str = "direct_occ") -> Path:
    case_directory = case_directory.expanduser().resolve()
    input_dir, results_dir = case_directory / "input", case_directory / "results"
    build_dir = case_directory / ".build"
    input_dir.mkdir(parents=True, exist_ok=True)
    results_dir.mkdir(parents=True, exist_ok=True)
    if operating_mode not in {"benchmark", "main_inlets", "branch_inlet"}:
        raise ValueError(f"Unsupported T-junction operating mode: {operating_mode}")
    if meshing_route not in {"direct_occ", "stl"}:
        raise ValueError(f"Unsupported T-junction meshing route: {meshing_route}")
    executable = build(build_dir, compiler)
    gmsh_executable = shutil.which(gmsh) or gmsh
    benchmark_path = SOURCE.parents[2] / "automatic_tetrahedra" / "benchmark.json"
    benchmark = json.loads(benchmark_path.read_text(encoding="utf-8"))
    geometry = benchmark["geometry"]
    stl = input_dir / "pipe_tee_oecd_nea_vattenfall.stl"
    subprocess.run(
        [str(executable), "--output", str(stl), "--gmsh", gmsh_executable,
         "--python", sys.executable, "--extractor", str(EXTRACTOR),
         "--mesh-size", str(mesh_size),
         "--main-diameter", str(geometry["main_pipe"]["inner_diameter"]),
         "--branch-diameter", str(geometry["branch_pipe"]["inner_diameter"]),
         "--main-inlet-length", str(geometry["main_pipe"]["inlet_length"]),
         "--main-outlet-length", str(geometry["main_pipe"]["outlet_length"]),
         "--branch-inlet-length", str(geometry["branch_pipe"]["inlet_length"]),
         "--operating-mode", operating_mode],
        cwd=ROOT, check=True,
    )
    mesh = results_dir / "pipe_tee_oecd_nea_vattenfall.msh"
    if meshing_route == "stl":
        from grid.stl_to_gmsh.generate_volume_mesh import generate as generate_volume_mesh

        report = generate_volume_mesh(
            stl, mesh, gmsh_executable, mesh_size, mesh_size
        )
        if report.get("validated") is not True:
            raise RuntimeError("Generic STL-derived T-junction volume meshing failed.")
    else:
        mesh_occ_volume(stl.with_suffix(".geo"), mesh, gmsh, geometry, operating_mode)
    return mesh


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case-directory", type=Path, required=True)
    parser.add_argument("--gmsh", required=True)
    parser.add_argument("--mesh-size", type=float, default=0.03)
    parser.add_argument("--compiler")
    parser.add_argument("--operating-mode", choices=("benchmark", "main_inlets", "branch_inlet"), default="benchmark")
    parser.add_argument("--meshing-route", choices=("direct_occ", "stl"), default="direct_occ")
    args = parser.parse_args()
    try:
        mesh = run_case(
            args.case_directory, gmsh=args.gmsh, mesh_size=args.mesh_size,
            compiler=args.compiler, operating_mode=args.operating_mode,
            meshing_route=args.meshing_route,
        )
    except (OSError, RuntimeError, ValueError, subprocess.CalledProcessError) as error:
        parser.exit(2, f"Native C++ T-junction Case failed: {error}\n")
    print(f"PASS: native C++ T-junction Case generated {mesh}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
