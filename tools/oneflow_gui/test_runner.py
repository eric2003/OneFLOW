"""Build commands for the existing OneFLOW unit and regression test runners."""

from __future__ import annotations

import os
import shlex
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

from .task_types import ExecutionLocation, TaskType, task_definition


@dataclass(frozen=True, slots=True)
class TestCommand:
    """One existing test command and the directory in which it must run."""

    label: str
    arguments: tuple[str, ...]
    working_directory: Path
    task_type: TaskType = TaskType.UNIT_TESTS


_PYTHON_UNIT_TESTS = (
    ("Workbench task type tests", ("-m", "unittest", "tools.oneflow_gui.test_task_types")),
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
    include_ctest: bool | None = None,
) -> list[TestCommand]:
    """Return the same CTest and Python test commands used by repository CI."""
    if build_config not in {"Debug", "Release"}:
        raise ValueError(f"Unsupported CMake build configuration: {build_config!r}")
    task = task_definition(TaskType.UNIT_TESTS)
    if task.execution_location is not ExecutionLocation.SOURCE_BUILD:
        raise RuntimeError("Unit tests must execute from the source build tree.")
    root = Path(repository_root).resolve()
    python_executable = python or sys.executable
    commands: list[TestCommand] = []
    build_directory = root / "build"
    ctest_is_available = (build_directory / "CTestTestfile.cmake").is_file()
    should_include_ctest = (
        ctest_is_available if include_ctest is None else include_ctest
    )
    if should_include_ctest:
        arguments = ["ctest", "--output-on-failure", "-C", build_config]
        commands.append(TestCommand("C++ unit tests (CTest)", tuple(arguments), build_directory))

    for label, arguments in _PYTHON_UNIT_TESTS:
        commands.append(
            TestCommand(label, (python_executable, *arguments), root)
        )
    return commands



def build_cmake_command(
    repository_root: Path,
    build_config: str = "Debug",
    *,
    rebuild: bool = False,
) -> TestCommand:
    """Build the selected single- or multi-configuration CMake build tree."""
    if build_config not in {"Debug", "Release"}:
        raise ValueError(f"Unsupported CMake build configuration: {build_config!r}")
    root = Path(repository_root).resolve()
    build_directory = root / "build"
    arguments = ["cmake", "--build", str(build_directory), "--config", build_config]
    if rebuild:
        arguments.append("--clean-first")
    label = f"CMake {'rebuild' if rebuild else 'build'} ({build_config})"
    return TestCommand(label, tuple(arguments), root)


def build_cmake_configure_command(
    repository_root: Path,
    build_config: str = "Debug",
) -> TestCommand:
    """Configure the repository build tree using the selected build configuration."""
    if build_config not in {"Debug", "Release"}:
        raise ValueError(f"Unsupported CMake build configuration: {build_config!r}")
    root = Path(repository_root).resolve()
    build_directory = root / "build"
    arguments = (
        "cmake",
        "..",
        "-D",
        "CMAKE_INSTALL_PREFIX=C:/dev/OneFLOW",
        "-D",
        "CMAKE_BUILD_TYPE=" + build_config,
    )
    return TestCommand(f"CMake configure ({build_config})", arguments, build_directory)


def cmake_build_tree_is_configured(build_directory: Path) -> bool:
    """Check for both CMake cache metadata and a generated build-system file."""
    build = Path(build_directory)
    has_cache = (build / "CMakeCache.txt").is_file()
    generated_files = ("build.ninja", "Makefile", "OneFLOW.sln")
    has_build_system = any((build / name).is_file() for name in generated_files)
    has_ctest_configuration = (build / "CTestTestfile.cmake").is_file()
    return has_cache and has_build_system and has_ctest_configuration


def build_default_test_commands(
    repository_root: Path,
    python: str | None = None,
    *,
    build_config: str = "Debug",
) -> list[TestCommand]:
    """Build the selected executable when needed, then run CTest and Python tests."""
    if build_config not in {"Debug", "Release"}:
        raise ValueError(f"Unsupported CMake build configuration: {build_config!r}")
    root = Path(repository_root).resolve()
    executable = root / "build" / "bin" / build_config / "OneFLOW.exe"
    commands: list[TestCommand] = []
    if not executable.is_file():
        if not cmake_build_tree_is_configured(root / "build"):
            commands.append(build_cmake_configure_command(root, build_config))
        commands.append(build_cmake_command(root, build_config))
        commands.extend(
            build_unit_test_commands(
                root,
                python,
                build_config=build_config,
                include_ctest=True,
            )
        )
        return commands
    return build_unit_test_commands(root, python, build_config=build_config)


def available_suites(repository_root: Path) -> list[str]:
    """List existing suite manifests without introducing a second suite format."""
    suite_directory = Path(repository_root).resolve() / "test" / "suites"
    if not suite_directory.is_dir():
        return []
    return sorted(
        path.name for path in suite_directory.glob("*.txt") if path.is_file()
    )



def stage_regression_tests(
    repository_root: Path,
    destination: Path,
) -> Path:
    """Replace the dedicated user workspace with a fresh copy of the test tree."""
    task = task_definition(TaskType.CFD_GRID_REGRESSION)
    if task.execution_location is not ExecutionLocation.USER_WORKSPACE:
        raise RuntimeError("CFD/grid regression must execute from a user workspace.")
    root = Path(repository_root).resolve()
    source_test_directory = root / "test"
    if not source_test_directory.is_dir():
        raise FileNotFoundError(f"Test directory not found: {source_test_directory}")

    target = Path(destination).expanduser().resolve()
    try:
        target.relative_to(root)
    except ValueError:
        pass
    else:
        raise ValueError("CFD regression workspace must be outside the source repository")

    if target == root or target in root.parents:
        raise ValueError("CFD regression workspace cannot contain the source repository")

    if target.exists():
        if target.is_dir() and not target.is_symlink():
            shutil.rmtree(target)
        else:
            target.unlink()
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copytree(source_test_directory, target)
    return target


def build_regression_command(
    repository_root: Path,
    executable: str | Path,
    suite_name: str,
    *,
    python: str | None = None,
    platform_name: str | None = None,
    test_directory: Path | None = None,
    validate_executable: bool = True,
) -> TestCommand:
    """Invoke test/test.py with an existing suite manifest and the platform MPI command."""
    root = Path(repository_root).resolve()
    regression_directory = (test_directory or (root / "test")).expanduser().resolve()
    suite_directory = (regression_directory / "suites").resolve()
    suite_path = (suite_directory / suite_name).resolve()
    if suite_path.parent != suite_directory or not suite_path.is_file():
        raise ValueError(f"Unknown regression suite: {suite_name!r}")
    executable_path = Path(executable).expanduser().resolve()
    if validate_executable and not executable_path.is_file():
        raise FileNotFoundError(f"OneFLOW executable not found: {executable_path}")

    platform_value = (platform_name or os.name).lower()
    is_windows = platform_value in {"nt", "windows"}
    mpi_command = "mpiexec -n 1" if is_windows else "mpirun -np 1"
    executable_argument = (
        subprocess.list2cmdline([str(executable_path)])
        if is_windows
        else shlex.quote(str(executable_path))
    )
    return TestCommand(
        f"CFD regression suite: {suite_name}",
        (
            python or sys.executable,
            "test.py",
            mpi_command,
            executable_argument,
            f"suites/{suite_name}",
        ),
        regression_directory,
        task_type=TaskType.CFD_GRID_REGRESSION,
    )
