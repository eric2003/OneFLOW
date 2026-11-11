#!/usr/bin/env python3
"""Stage and execute the registered STL-to-Gmsh development Cases."""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

REPOSITORY_ROOT = Path(__file__).resolve().parents[3]
NATIVE_CASES = {
    "square/single_quad/cpp": ("square/single_quad/cpp", "square.msh"),
    "square/two_triangles/cpp": (
        "square/two_triangles/cpp", "square_two_triangles.msh"
    ),
    "square/automatic_triangles/cpp": (
        "square/automatic_triangles/cpp", "square_automatic_triangles.msh"
    ),
    "square_2x2/cpp": ("square_2x2/cpp", "square_2x2.msh"),
    "two_square_cells/cpp": ("two_square_cells/cpp", "two_square_cells.msh"),
    "circle_external_flow/cpp": ("circle_external_flow/cpp", "circle_external_flow.msh"),
    "cube/single_hex/cpp": ("cube/single_hex/cpp", "cube_single_hex.msh"),
    "cube/six_tetrahedra/cpp": ("cube/six_tetrahedra/cpp", "cube_six_tetrahedra.msh"),
    "cube/automatic_tetrahedra/cpp": ("cube/automatic_tetrahedra/cpp", "cube_automatic_tetrahedra.msh"),
}
PYTHON_CASES = {
    "square/single_quad/python": "square/single_quad/python",
    "square/two_triangles/python": "square/two_triangles/python",
    "square/automatic_triangles/python": "square/automatic_triangles/python",
    "two_square_cells/python": "two_square_cells/python",
    "circle_external_flow/python": "circle_external_flow/python",
    "channel_internal_flow/python": "channel_internal_flow/python",
    "cube/single_hex/python": "cube/single_hex/python",
    "cube/six_tetrahedra/python": "cube/six_tetrahedra/python",
    "cube/automatic_tetrahedra/python": "cube/automatic_tetrahedra/python",
}


def find_oneflow_executable(explicit: str | None) -> Path:
    if explicit:
        path = Path(explicit).expanduser().resolve()
        if not path.is_file():
            raise FileNotFoundError(f"OneFLOW executable does not exist: {path}")
        return path
    candidates = (
        "build/bin/Debug/OneFLOW.exe",
        "build/bin/Release/OneFLOW.exe",
        "build/bin/OneFLOW.exe",
        "build/bin/Debug/OneFLOW",
        "build/bin/Release/OneFLOW",
        "build/bin/OneFLOW",
    )
    for candidate in candidates:
        path = REPOSITORY_ROOT / candidate
        if path.is_file():
            return path
    raise FileNotFoundError(
        "OneFLOW executable was not found; pass --oneflow-executable"
    )


def find_gmsh_executable(explicit: str) -> Path:
    candidate = Path(explicit).expanduser()
    if candidate.is_file():
        return candidate.resolve()
    discovered = shutil.which(explicit)
    if discovered:
        return Path(discovered).resolve()
    raise FileNotFoundError(f"Gmsh executable does not exist or is not on PATH: {explicit}")


def stage_case(source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists():
        shutil.rmtree(destination)
    shutil.copytree(source, destination)


def execute(command: list[str], case_directory: Path) -> None:
    result = subprocess.run(
        command,
        cwd=case_directory,
        text=True,
        encoding="utf-8",
        errors="replace",
        capture_output=True,
        check=False,
    )
    results = case_directory / "results"
    results.mkdir(parents=True, exist_ok=True)
    log = results / "development_case.log"
    log.write_text(
        "$ " + subprocess.list2cmdline(command) + "\n\n"
        + result.stdout + result.stderr,
        encoding="utf-8",
        newline="\n",
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"Command failed with exit code {result.returncode}; see {log}\n"
            f"{result.stdout}\n{result.stderr}"
        )


def run_all(stage_root: Path, oneflow: Path, gmsh: str) -> list[tuple[str, Path]]:
    outputs: list[tuple[str, Path]] = []
    for case_id, (directory, mesh_name) in NATIVE_CASES.items():
        source = REPOSITORY_ROOT / "grid/stl_to_gmsh" / directory
        destination = stage_root / case_id
        stage_case(source, destination)
        execute([str(oneflow), "r", str(destination)], destination)
        mesh = destination / "results" / mesh_name
        report = destination / "results" / f"{Path(mesh_name).stem}.stl-to-gmsh-report.json"
        if not mesh.is_file() or not report.is_file():
            raise RuntimeError(
                f"{case_id}: expected mesh/report outputs are missing under "
                f"{destination / 'results'}"
            )
        outputs.append((case_id, destination))

    python_runner = REPOSITORY_ROOT / "grid/stl_to_gmsh/run_python_case.py"
    for case_id, directory in PYTHON_CASES.items():
        source = REPOSITORY_ROOT / "grid/stl_to_gmsh" / directory
        destination = stage_root / case_id
        stage_case(source, destination)
        execute(
            [
                sys.executable,
                str(python_runner),
                "--case-directory",
                str(destination),
                "--gmsh",
                gmsh,
            ],
            destination,
        )
        expected = json.loads(
            (destination / "expected.json").read_text(encoding="utf-8")
        )
        mesh_name = expected["artifact"]["mesh_file"]
        report_name = f"{Path(mesh_name).stem}.stl-to-gmsh-report.json"
        mesh = destination / "results" / mesh_name
        report = destination / "results" / report_name
        if not mesh.is_file() or not report.is_file():
            raise RuntimeError(
                f"{case_id}: expected mesh/report outputs are missing under "
                f"{destination / 'results'}"
            )
        outputs.append((case_id, destination))
    return outputs


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--oneflow-executable",
        help="OneFLOW executable; defaults to a build/bin Debug/Release candidate",
    )
    parser.add_argument(
        "--gmsh",
        default=shutil.which("gmsh"),
        help="Gmsh executable (default: search PATH)",
    )
    parser.add_argument(
        "--stage-root",
        type=Path,
        default=Path.home() / ".oneflow/mesh-layout-test/grid/stl_to_gmsh",
        help="staging directory outside the source tree",
    )
    args = parser.parse_args()
    try:
        if not args.gmsh:
            raise FileNotFoundError("Gmsh executable was not found; pass --gmsh")
        oneflow = find_oneflow_executable(args.oneflow_executable)
        gmsh_path = find_gmsh_executable(args.gmsh)
        stage_root = args.stage_root.expanduser().resolve()
        os.environ["PATH"] = (
            str(gmsh_path.parent) + os.pathsep + os.environ.get("PATH", "")
        )
        try:
            stage_root.relative_to(REPOSITORY_ROOT)
        except ValueError:
            pass
        else:
            raise ValueError(
                "The staging directory must be outside the OneFLOW source repository."
            )

        outputs = run_all(stage_root, oneflow, str(gmsh_path))
    except (OSError, RuntimeError, ValueError) as error:
        parser.exit(2, f"STL-to-Gmsh development run failed: {error}\n")

    if outputs:
        print("STL-to-Gmsh development Cases passed:")
        for case_id, directory in outputs:
            print(f"  {case_id}: {directory / 'results'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
