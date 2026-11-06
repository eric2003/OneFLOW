"""Manual GUI runner for the experimental Gmsh-to-CGNS bootstrap Case."""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


def resolve_executable(value: str) -> Path:
    candidate = Path(value).expanduser()
    if candidate.is_file():
        return candidate.resolve()
    located = shutil.which(value)
    if located:
        return Path(located).resolve()
    raise FileNotFoundError(f"OneFLOW executable not found: {value}")


def run_case(executable: Path, repository_root: Path) -> int:
    source_case = repository_root / "grid" / "gmsh_to_cgns"
    if not source_case.is_dir():
        raise FileNotFoundError(f"Gmsh Case not found: {source_case}")

    with tempfile.TemporaryDirectory(prefix="oneflow_gmsh_to_cgns_") as workspace:
        case_dir = Path(workspace) / "gmsh_to_cgns"
        shutil.copytree(source_case, case_dir)
        command = [str(executable), "r", str(case_dir)]
        print("Running experimental Gmsh Case:", " ".join(command), flush=True)
        completed = subprocess.run(command, cwd=repository_root, check=False)
        if completed.returncode != 0:
            print(
                f"FAIL: OneFLOW exited with code {completed.returncode}",
                file=sys.stderr,
            )
            return completed.returncode or 1

        artifact = case_dir / "results" / "square.cgns"
        if not artifact.is_file():
            print(f"FAIL: expected output is missing: {artifact}", file=sys.stderr)
            return 1

        content = artifact.read_text(encoding="utf-8")
        required_markers = (
            "ONEFLOW_CGNS_PLACEHOLDER 1",
            "status=NOT_A_CGNS_FILE",
            "source_format=gmsh_msh_2_ascii",
            "node_count=4",
            "element_count=2",
        )
        missing = [marker for marker in required_markers if marker not in content]
        if missing:
            print(
                "FAIL: output does not contain expected conversion evidence: "
                + ", ".join(missing),
                file=sys.stderr,
            )
            return 1

        print("PASS: Gmsh Case executed and produced the explicitly marked placeholder.")
        print("NOTE: this artifact is not a standards-compliant CGNS file.")
        print("Workspace cleanup: temporary Case and outputs will be removed.")
        return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repository-root", required=True)
    parser.add_argument("--executable", required=True)
    args = parser.parse_args()
    try:
        root = Path(args.repository_root).expanduser().resolve()
        executable = resolve_executable(args.executable)
        return run_case(executable, root)
    except (OSError, UnicodeError) as error:
        print(f"Experimental Gmsh Case failed: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
