#!/usr/bin/env python3
"""Compile and run the C++17 T-junction geometry driver, then use shared Gmsh meshing."""
from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[6]
SOURCE = ROOT / "grid/stl_to_gmsh/3d/pipe_tee/oecd_nea_vattenfall/geometry/cpp/generate_pipe_tee.cpp"
EXTRACTOR = ROOT / "grid/stl_to_gmsh/3d/pipe_tee/oecd_nea_vattenfall/geometry/python/generate_pipe_tee.py"
VOLUME_MESHER = ROOT / "grid/stl_to_gmsh/generate_volume_mesh.py"


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
             compiler: str | None = None) -> Path:
    case_directory = case_directory.expanduser().resolve()
    input_dir, results_dir = case_directory / "input", case_directory / "results"
    build_dir = case_directory / ".build"
    input_dir.mkdir(parents=True, exist_ok=True)
    results_dir.mkdir(parents=True, exist_ok=True)
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
         "--branch-inlet-length", str(geometry["branch_pipe"]["inlet_length"])],
        cwd=ROOT, check=True,
    )
    mesh = results_dir / "pipe_tee_oecd_nea_vattenfall.msh"
    subprocess.run(
        [sys.executable, str(VOLUME_MESHER), "--stl", str(stl), "--output", str(mesh),
         "--gmsh", gmsh, "--mesh-size", str(mesh_size),
         "--boundary-name", "wall", "--volume-name", "fluid_domain"],
        cwd=ROOT, check=True,
    )
    return mesh


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case-directory", type=Path, required=True)
    parser.add_argument("--gmsh", required=True)
    parser.add_argument("--mesh-size", type=float, default=0.03)
    parser.add_argument("--compiler")
    args = parser.parse_args()
    try:
        mesh = run_case(args.case_directory, gmsh=args.gmsh, mesh_size=args.mesh_size, compiler=args.compiler)
    except (OSError, RuntimeError, subprocess.CalledProcessError) as error:
        parser.exit(2, f"Native C++ T-junction Case failed: {error}\n")
    print(f"PASS: native C++ T-junction Case generated {mesh}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
