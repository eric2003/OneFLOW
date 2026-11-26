#!/usr/bin/env python3
"""Run the native C++ pipe-bend geometry producer through the shared Gmsh mesher."""
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

REPOSITORY_ROOT = Path(__file__).resolve().parents[6]
GEOMETRY_SOURCE = (
    REPOSITORY_ROOT
    / "grid/stl_to_gmsh/3d/pipe_bend/90_degree/geometry/cpp/generate_pipe_bend.cpp"
)
VOLUME_MESHER = REPOSITORY_ROOT / "grid/stl_to_gmsh/generate_volume_mesh.py"


def discover_cpp_compiler() -> str:
    """Find a C++17 compiler, honoring the user's explicit CXX selection."""
    requested = os.environ.get("CXX", "").strip()
    if requested:
        compiler = shutil.which(requested)
        if compiler:
            return compiler
        raise RuntimeError(f"CXX compiler was configured but not found: {requested}")
    for name in ("c++", "g++", "clang++", "cl"):
        compiler = shutil.which(name)
        if compiler:
            return compiler
    raise RuntimeError(
        "A C++17 compiler is required for the native pipe-bend implementation. "
        "Install GCC/Clang, or run from a Visual Studio Developer Command Prompt."
    )


def build_cpp_generator(
    build_directory: Path, *, compiler: str | None = None
) -> Path:
    """Compile the checked-in generator into the external Case workspace."""
    compiler = compiler or discover_cpp_compiler()
    build_directory = build_directory.expanduser().resolve()
    build_directory.mkdir(parents=True, exist_ok=True)
    executable = build_directory / (
        "pipe_bend_generator.exe" if os.name == "nt" else "pipe_bend_generator"
    )
    if Path(compiler).name.lower() in {"cl", "cl.exe"}:
        command = [
            compiler, "/nologo", "/std:c++17", "/EHsc", "/O2",
            str(GEOMETRY_SOURCE), f"/Fe:{executable}",
        ]
    else:
        command = [
            compiler, "-std=c++17", "-O2", str(GEOMETRY_SOURCE),
            "-o", str(executable),
        ]
    print("Building native geometry:", subprocess.list2cmdline(command), flush=True)
    subprocess.run(command, cwd=build_directory, check=True)
    return executable


def run_case(
    case_directory: Path,
    *,
    gmsh: str,
    mesh_size: float = 0.20,
    compiler: str | None = None,
) -> Path:
    """Generate native C++ geometry, then use the same generic volume mesher."""
    case_directory = case_directory.expanduser().resolve()
    input_directory = case_directory / "input"
    results_directory = case_directory / "results"
    build_directory = case_directory / ".build"
    input_directory.mkdir(parents=True, exist_ok=True)
    results_directory.mkdir(parents=True, exist_ok=True)

    executable = build_cpp_generator(build_directory, compiler=compiler)
    stl = input_directory / "pipe_bend_90_degree.stl"
    mesh = results_directory / "pipe_bend_90_degree.msh"
    geometry_command = [str(executable), "--output", str(stl)]
    print("Running native geometry:", subprocess.list2cmdline(geometry_command), flush=True)
    subprocess.run(geometry_command, cwd=REPOSITORY_ROOT, check=True)

    mesh_command = [
        sys.executable, str(VOLUME_MESHER),
        "--stl", str(stl), "--output", str(mesh),
        "--gmsh", gmsh, "--mesh-size", str(mesh_size),
        "--boundary-name", "wall", "--volume-name", "fluid_domain",
    ]
    print("Running shared volume mesher:", subprocess.list2cmdline(mesh_command), flush=True)
    subprocess.run(mesh_command, cwd=REPOSITORY_ROOT, check=True)
    return mesh


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case-directory", type=Path, required=True)
    parser.add_argument("--gmsh", required=True)
    parser.add_argument("--mesh-size", type=float, default=0.20)
    parser.add_argument("--compiler", help="C++17 compiler; defaults to CXX or PATH discovery")
    args = parser.parse_args()
    try:
        mesh = run_case(
            args.case_directory, gmsh=args.gmsh,
            mesh_size=args.mesh_size, compiler=args.compiler,
        )
    except (OSError, RuntimeError, subprocess.CalledProcessError) as error:
        parser.exit(2, f"Native pipe-bend Case failed: {error}\n")
    print(f"PASS: native C++ pipe-bend Case generated {mesh}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
