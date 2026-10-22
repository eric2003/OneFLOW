"""Build commands for the existing OneFLOW unit and regression test runners."""

from __future__ import annotations

import os
import shlex
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True, slots=True)
class TestCommand:
    """One existing test command and the directory in which it must run."""

    label: str
    arguments: tuple[str, ...]
    working_directory: Path


_PYTHON_UNIT_TESTS = (
    ("Mesh workflow tests", ("-m", "unittest", "discover", "-s", "tools", "-p", "test_mesh_workflow.py")),
    ("Mesh export tests", ("-m", "unittest", "discover", "-s", "tools/mesh_export", "-p", "test_*.py")),
    ("Mesh case catalog tests", ("tools/test_mesh_cases.py",)),
    ("Flat-plate layout compiler tests", ("-m", "unittest", "tests.tools.test_compile_flat_plate")),
    ("Explicit layout topology tests", ("-m", "unittest", "tests.tools.test_layout_topology")),
    ("Desktop service tests", ("-m", "unittest", "tools.oneflow_gui.test_services")),
    ("Project repository tests", ("-m", "unittest", "tools.oneflow_gui.test_project_repository")),
    ("User data path tests", ("-m", "unittest", "tools.oneflow_gui.test_user_data_paths")),
    ("Workflow persistence tests", ("-m", "unittest", "tools.oneflow_gui.test_workflow_repository")),
    ("Workflow command tests", ("-m", "unittest", "tools.oneflow_gui.test_workflow_execution")),
    ("Workflow dependency tests", ("-m", "unittest", "tools.oneflow_gui.test_workflow_dependencies")),
    ("Workflow plan tests", ("-m", "unittest", "tools.oneflow_gui.test_workflow_plan")),
    ("Workflow graph runner tests", ("-m", "unittest", "tools.oneflow_gui.test_workflow_graph_runner")),
    ("Workflow process tests", ("-m", "unittest", "tools.oneflow_gui.test_workflow_process")),
    ("Runtime override tests", ("-m", "unittest", "tools.oneflow_gui.test_workflow_runtime_overrides")),
    ("Workflow case staging tests", ("-m", "unittest", "tools.oneflow_gui.test_workflow_case_runtime")),
    ("GUI test command builder tests", ("-m", "unittest", "tools.oneflow_gui.test_test_runner")),
)


def build_unit_test_commands(
    repository_root: Path,
    python: str | None = None,
    *,
    build_config: str = "Debug",
) -> list[TestCommand]:
    """Return the same CTest and Python test commands used by repository CI."""
    if build_config not in {"Debug", "Release"}:
        raise ValueError(f"Unsupported CMake build configuration: {build_config!r}")
    root = Path(repository_root).resolve()
    python_executable = python or sys.executable
    commands: list[TestCommand] = []
    build_directory = root / "build"
    if (build_directory / "CTestTestfile.cmake").is_file():
        arguments = ["ctest", "--output-on-failure", "-C", build_config]
        commands.append(TestCommand("C++ unit tests (CTest)", tuple(arguments), build_directory))

    for label, arguments in _PYTHON_UNIT_TESTS:
        commands.append(
            TestCommand(label, (python_executable, *arguments), root)
        )
    return commands


def available_suites(repository_root: Path) -> list[str]:
    """List existing suite manifests without introducing a second suite format."""
    suite_directory = Path(repository_root).resolve() / "test" / "suites"
    if not suite_directory.is_dir():
        return []
    return sorted(
        path.name for path in suite_directory.glob("*.txt") if path.is_file()
    )


def build_regression_command(
    repository_root: Path,
    executable: str | Path,
    suite_name: str,
    *,
    python: str | None = None,
    platform_name: str | None = None,
) -> TestCommand:
    """Invoke test/test.py with an existing suite manifest and the platform MPI command."""
    root = Path(repository_root).resolve()
    suite_directory = (root / "test" / "suites").resolve()
    suite_path = (suite_directory / suite_name).resolve()
    if suite_path.parent != suite_directory or not suite_path.is_file():
        raise ValueError(f"Unknown regression suite: {suite_name!r}")
    executable_path = Path(executable).expanduser().resolve()
    if not executable_path.is_file():
        raise FileNotFoundError(f"OneFLOW executable not found: {executable_path}")

    platform_value = (platform_name or os.name).lower()
    is_windows = platform_value in {"nt", "windows"}
    mpi_command = "mpiexec -n 1" if is_windows else "mpirun -np 1"
    executable_argument = (
        subprocess.list2cmdline([str(executable_path)])
        if is_windows
        else shlex.quote(str(executable_path))
    )
    test_directory = root / "test"
    return TestCommand(
        f"CFD regression suite: {suite_name}",
        (
            python or sys.executable,
            "test.py",
            mpi_command,
            executable_argument,
            f"suites/{suite_name}",
        ),
        test_directory,
    )
