#!/usr/bin/env python3
"""Run the native C++ T-junction producer through the generic STL mesher."""
from __future__ import annotations

import argparse
from pathlib import Path

from run_pipe_tee_cpp_case import run_case


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case-directory", type=Path, required=True)
    parser.add_argument("--gmsh", required=True)
    parser.add_argument("--mesh-size", type=float, default=0.03)
    parser.add_argument("--compiler")
    parser.add_argument(
        "--operating-mode",
        choices=("benchmark", "main_inlets", "branch_inlet"),
        default="main_inlets",
    )
    args = parser.parse_args()
    try:
        mesh = run_case(
            args.case_directory,
            gmsh=args.gmsh,
            mesh_size=args.mesh_size,
            compiler=args.compiler,
            operating_mode=args.operating_mode,
            meshing_route="stl",
        )
    except (OSError, RuntimeError, ValueError) as error:
        parser.exit(2, f"Native C++ T-junction STL Case failed: {error}\\n")
    print(f"PASS: native C++ STL-derived T-junction mesh generated {mesh}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
