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


def format_command_arguments(arguments: tuple[str, ...] | list[str], *, platform_name: str | None = None) -> str:
    """Render an argument vector for logs using the host command-line syntax."""
    platform_value = (platform_name or os.name).lower()
    if platform_value in {"nt", "windows"}:
        return subprocess.list2cmdline(list(arguments))
    return shlex.join(str(argument) for argument in arguments)


def oneflow_executable_path(
    repository_root: Path,
    build_config: str = "Debug",
    *,
    platform_name: str | None = None,
) -> Path:
    """Return the conventional build output for the current platform."""
    if build_config not in {"Debug", "Release"}:
        raise ValueError(f"Unsupported CMake build configuration: {build_config!r}")
    root = Path(repository_root).resolve()
    platform_value = (platform_name or os.name).lower()
    build_directory = root / "build"
    if platform_value in {"nt", "windows"}:
        return build_directory / "bin" / build_config / "OneFLOW.exe"

    single_config_executable = build_directory / "bin" / "OneFLOW"
    if single_config_executable.is_file():
        return single_config_executable

    multi_config_cache = build_directory / "CMakeCache.txt"
    if multi_config_cache.is_file():
        cache_text = multi_config_cache.read_text(encoding="utf-8", errors="replace")
        if "CMAKE_CONFIGURATION_TYPES:" in cache_text:
            return build_directory / "bin" / build_config / "OneFLOW"
    return single_config_executable


def summarize_test_output(
    output_lines: list[str], *, include_failed_cases: bool = True
) -> tuple[str, ...]:
    """Extract framework test counts and actionable failure details from output."""
    import re

    output = "\n".join(output_lines)
    summary: list[str] = []
    unittest_results = re.findall(
        r"Ran (\d+) tests?[^\n]*(?:\n[ \t]*)?(?:\n[ \t]*)?"
        r"(FAILED \(([^)]*)\)|OK(?: \(([^)]*)\))?)",
        output,
    )
    ctest_matches = re.findall(
        r"(\d+)% tests passed, (\d+) tests failed out of (\d+)",
        output,
    )
    if unittest_results:
        total = sum(int(item[0]) for item in unittest_results)
        failed = 0
        skipped = 0
        for _, _, failed_status, ok_status in unittest_results:
            status = failed_status or ok_status or ""

            def status_count(name: str) -> int:
                match = re.search(rf"{name}=(\d+)", status)
                return int(match.group(1)) if match else 0

            failed += status_count("failures") + status_count("errors")
            skipped += status_count("skipped")
        passed = max(total - failed - skipped, 0)
        summary.append(
            f"Python unittest cases: {total} run, {passed} passed, "
            f"{failed} failed, {skipped} skipped."
        )
    if ctest_matches:
        total = sum(int(item[2]) for item in ctest_matches)
        failed = sum(int(item[1]) for item in ctest_matches)
        skipped_block = re.search(
            r"The following tests did not run:\s*((?:\n[ \t]+\d+\s+-[^\n]*)+)",
            output,
        )
        skipped = (
            len(re.findall(r"^\s*\d+\s+-", skipped_block.group(1), re.MULTILINE))
            if skipped_block
            else 0
        )
        passed = max(total - failed - skipped, 0)
        summary.append(
            f"CTest cases: {total} run, {passed} passed, "
            f"{failed} failed, {skipped} skipped."
        )

    failure_markers = (
        "FAIL:",
        "ERROR:",
        "FAILED:",
        "The following tests FAILED:",
        "Traceback (most recent call last):",
        "fatal error",
        "error ",
        "AssertionError:",
    )
    failure_start = next(
        (
            index
            for index, line in enumerate(output_lines)
            if line.lstrip().startswith(failure_markers)
        ),
        None,
    )
    failed_unittest_cases = []
    for line in output_lines:
        match = re.match(r"^\s*(FAIL|ERROR):\s+(.+?)\s*$", line)
        if match:
            case_name = match.group(2)
            if case_name not in failed_unittest_cases:
                failed_unittest_cases.append(case_name)
    if failed_unittest_cases and include_failed_cases:
        summary.append("Failed unittest cases:")
        summary.extend(f"  - {case_name}" for case_name in failed_unittest_cases)

    if failure_start is not None:
        details = output_lines[failure_start:]
        summary.append("Failure details:")
        summary.extend(details[:120])
        if len(details) > 120:
            summary.append(f"... {len(details) - 120} additional lines omitted.")
    return tuple(summary)


def format_final_test_status(passed: bool, report: str) -> str:
    """Format command and framework case counts on the final status line."""
    import re

    status = "PASS" if passed else "FAIL"
    details: list[str] = []
    command_counts = re.search(
        r"Test command summary:\s*\n\s*Total: (\d+)\s*\n"
        r"\s*Passed: (\d+)\s*\n\s*Failed: (\d+)\s*\n"
        r"\s*Not run: (\d+)",
        report,
    )
    if command_counts:
        total, command_passed, command_failed, not_run = map(
            int, command_counts.groups()
        )
        details.append(
            f"Commands: {command_passed} passed, {command_failed} failed, "
            f"{not_run} not run ({total} total)"
        )

    case_patterns = (
        (
            "Python unittest",
            r"Python unittest cases: (\d+) run, (\d+) passed, "
            r"(\d+) failed, (\d+) skipped\.",
        ),
        (
            "CTest",
            r"CTest cases: (\d+) run, (\d+) passed, "
            r"(\d+) failed, (\d+) skipped\.",
        ),
    )
    for label, pattern in case_patterns:
        matches = re.findall(pattern, report)
        if not matches:
            continue
        counts = [sum(int(item[index]) for item in matches) for index in range(4)]
        total, case_passed, case_failed, skipped = counts
        details.append(
            f"{label}: {case_passed} passed, {case_failed} failed, "
            f"{skipped} skipped ({total} total)"
        )

    suffix = " | " + " | ".join(details) if details else ""
    return f"FINAL STATUS: {status}{suffix}"


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
    ("Gmsh to OneFLOW case tests", ("tools/test_gmsh_to_oneflow_case.py",)),
    ("Planar STL to Gmsh tests", ("-m", "unittest", "tools.stl_to_gmsh.test_stl_to_gmsh")),
    ("MeshSnapshot contract tests", ("-m", "unittest", "tools.mesh_contracts.test_snapshot")),
    ("Unstructured CGNS writer architecture tests", ("-m", "unittest", "tools.mesh_contracts.test_unstructured_cgns_writer_architecture")),
    ("Flat-plate layout compiler tests", ("-m", "unittest", "tests.tools.test_compile_flat_plate")),
    ("Explicit layout topology tests", ("-m", "unittest", "tests.tools.test_layout_topology")),
    ("Desktop service tests", ("-m", "unittest", "tools.oneflow_gui.test_services")),
    ("Project repository tests", ("-m", "unittest", "tools.oneflow_gui.test_project_repository")),
    ("Case script discovery tests", ("-m", "unittest", "tools.oneflow_gui.test_script_file_discovery")),
    ("Project window tests", ("-m", "unittest", "tools.oneflow_gui.test_project_window")),
    ("Parameter document tests", ("-m", "unittest", "tools.oneflow_gui.test_parameter_document")),
    ("Logical parameter model tests", ("-m", "unittest", "tools.oneflow_gui.test_parameter_model")),
    ("File-independent parameter resolution tests", ("-m", "unittest", "tools.oneflow_gui.test_parameter_resolution")),
    ("Declarative Case profile tests", ("-m", "unittest", "tools.oneflow_gui.test_case_profile")),
    ("Runtime parameter consumer evidence tests", ("-m", "unittest", "tools.oneflow_gui.test_parameter_consumer_catalog")),
    ("CFD nomenclature knowledge tests", ("-m", "unittest", "tools.verification.test_nomenclature_catalog")),
    ("Parameter catalog database tests", ("-m", "unittest", "tools.oneflow_gui.test_parameter_catalog_store")),
    ("User data path tests", ("-m", "unittest", "tools.oneflow_gui.test_user_data_paths")),
    ("Application settings tests", ("-m", "unittest", "tools.oneflow_gui.test_application_settings")),
    ("Workflow persistence tests", ("-m", "unittest", "tools.oneflow_gui.test_workflow_repository")),
    ("Workspace persistence tests", ("-m", "unittest", "tools.oneflow_gui.test_workspace_repository")),
    ("Workspace execution adapter tests", ("-m", "unittest", "tools.oneflow_gui.test_workspace_workflow_adapter")),
    ("Workflow command tests", ("-m", "unittest", "tools.oneflow_gui.test_workflow_execution")),
    ("Workflow dependency tests", ("-m", "unittest", "tools.oneflow_gui.test_workflow_dependencies")),
    ("Workflow plan tests", ("-m", "unittest", "tools.oneflow_gui.test_workflow_plan")),
    ("Workflow graph runner tests", ("-m", "unittest", "tools.oneflow_gui.test_workflow_graph_runner")),
    ("Workflow process tests", ("-m", "unittest", "tools.oneflow_gui.test_workflow_process")),
    ("Runtime override tests", ("-m", "unittest", "tools.oneflow_gui.test_workflow_runtime_overrides")),
    ("Workflow case staging tests", ("-m", "unittest", "tools.oneflow_gui.test_workflow_case_runtime")),
    ("GUI test command builder tests", ("-m", "unittest", "tools.oneflow_gui.test_test_runner")),
    ("Structured test diagnostic tests", ("-m", "unittest", "tools.oneflow_gui.test_test_diagnostics")),
    ("Numerical verification catalog tests", ("-m", "unittest", "tools.verification.test_curve_catalog")),
    ("Vencat analytic reference tests", ("-m", "unittest", "tools.verification.test_vencat_reference")),
    ("Combined contract catalog tests", ("-m", "unittest", "tools.verification.test_contract_catalog")),
    ("Verification execution record tests", ("-m", "unittest", "tools.verification.test_execution_records")),
    ("Python execution capture tests", ("-m", "unittest", "tools.oneflow_gui.test_execution_capture")),
    ("Workflow execution evidence tests", ("-m", "unittest", "tools.oneflow_gui.test_workflow_execution_capture")),
    ("Case execution evidence tests", ("-m", "unittest", "tools.oneflow_gui.test_case_execution_capture")),
    ("Case contract residual runner tests", ("-m", "unittest", "tools.oneflow_gui.test_case_contract_run")),
    ("Test worker execution evidence tests", ("-m", "unittest", "tools.oneflow_gui.test_test_dialog_execution_capture")),
    ("CMake test catalog tests", ("-m", "unittest", "tools.verification.test_cmake_test_catalog")),
    ("Configured CTest discovery tests", ("-m", "unittest", "tools.verification.test_ctest_discovery")),
    ("Contract Workbench dialog tests", ("-m", "unittest", "tools.oneflow_gui.test_contract_workbench_dialog")),
    ("Source editor process tests", ("-m", "unittest", "tools.oneflow_gui.test_source_editor")),
    ("Repository test inventory discovery tests", ("-m", "unittest", "tools.oneflow_gui.test_test_inventory")),
    ("Executable discovery tests", ("-m", "unittest", "tools.oneflow_gui.test_executable_discovery")),
    ("Process control tests", ("-m", "unittest", "tools.oneflow_gui.test_process_control")),
    ("Case runner process tests", ("-m", "unittest", "tools.oneflow_gui.test_runner_process")),
    ("Mesh layout editor document tests", ("-m", "unittest", "tests.tools.test_editor_document")),
    ("Rectangle layout compiler tests", ("-m", "unittest", "tests.tools.test_compile_rectangle")),
    ("Geometry analytic contract tests", ("-m", "unittest", "discover", "-s", "tools", "-p", "test_analytic_contracts.py")),
    ("Geometry scenario persistence tests", ("-m", "unittest", "discover", "-s", "tools", "-p", "test_scenario_contract_persistence.py")),
    ("Geometry scenario import tests", ("-m", "unittest", "discover", "-s", "tools", "-p", "test_scenario_import.py")),
    ("Geometry workbench fixture tests", ("-m", "unittest", "discover", "-s", "tools/geometry_workbench", "-p", "test_workbench_stub.py")),
    ("Background project runner tests", ("-m", "unittest", "tools.oneflow_gui.test_background_project_runner")),
    ("Background task host tests", ("-m", "unittest", "tools.oneflow_gui.test_background_task_host")),
    ("Background task manager tests", ("-m", "unittest", "tools.oneflow_gui.test_background_task_manager")),
    ("Background task registry tests", ("-m", "unittest", "tools.oneflow_gui.test_background_task_registry")),
    ("Geometry contract bridge tests", ("-m", "unittest", "tools.oneflow_gui.test_geometry_contract_bridge")),
    ("Mesh layout project tests", ("-m", "unittest", "tools.oneflow_gui.test_mesh_layout_projects")),
    ("Project inventory tests", ("-m", "unittest", "tools.oneflow_gui.test_project_inventory")),
    ("Project template tests", ("-m", "unittest", "tools.oneflow_gui.test_project_templates")),
    ("RAE2822 project template tests", ("-m", "unittest", "tools.oneflow_gui.test_rae2822_project_template")),
    ("SA closure curve tests", ("-m", "unittest", "discover", "-s", "tools/verification", "-p", "test_sa_closure_curves.py")),
    ("Tkinter runner tests", ("gui/python/tkinter/test_oneflow_runner.py",)),
)

_PYTHON_TEST_EXCLUSIONS = {
    "tools/oneflow_gui/test_main_window.py": "Qt widget/queue tests; keep in a dedicated headless GUI test group.",
    "tools/oneflow_gui/test_new_project_dialog.py": "Qt dialog tests; keep in a dedicated headless GUI test group.",
    "tools/oneflow_gui/test_workspace_dialog.py": "Qt dialog tests; keep in a dedicated headless GUI test group.",
    "tools/verification/test_latex_preview.py": "Optional Matplotlib rendering dependency; run in the visualization-enabled test environment.",
}


def build_ctest_test_commands(
    repository_root: Path,
    *,
    build_config: str = "Debug",
) -> list[TestCommand]:
    """Return the CTest command when the source build tree has test metadata."""
    if build_config not in {"Debug", "Release"}:
        raise ValueError(f"Unsupported CMake build configuration: {build_config!r}")
    task = task_definition(TaskType.UNIT_TESTS)
    if task.execution_location is not ExecutionLocation.SOURCE_BUILD:
        raise RuntimeError("Unit tests must execute from the source build tree.")
    root = Path(repository_root).resolve()
    build_directory = root / "build"
    if not (build_directory / "CTestTestfile.cmake").is_file():
        return []
    return [
        TestCommand(
            "C++ unit tests (CTest)",
            ("ctest", "--output-on-failure", "-C", build_config),
            build_directory,
        )
    ]


def build_python_test_commands(
    repository_root: Path,
    python: str | None = None,
) -> list[TestCommand]:
    """Return the explicit Python unittest command registry."""
    task = task_definition(TaskType.UNIT_TESTS)
    if task.execution_location is not ExecutionLocation.SOURCE_BUILD:
        raise RuntimeError("Unit tests must execute from the source build tree.")
    root = Path(repository_root).resolve()
    python_executable = python or sys.executable
    return [
        TestCommand(label, (python_executable, *arguments), root)
        for label, arguments in _PYTHON_UNIT_TESTS
    ]


def build_unit_test_commands(
    repository_root: Path,
    python: str | None = None,
    *,
    build_config: str = "Debug",
    include_ctest: bool | None = None,
) -> list[TestCommand]:
    """Return the combined CTest and Python test commands for legacy callers."""
    if build_config not in {"Debug", "Release"}:
        raise ValueError(f"Unsupported CMake build configuration: {build_config!r}")
    root = Path(repository_root).resolve()
    commands = build_ctest_test_commands(root, build_config=build_config)
    if include_ctest is False:
        commands = []
    elif include_ctest is True and not commands:
        build_directory = root / "build"
        commands = [
            TestCommand(
                "C++ unit tests (CTest)",
                ("ctest", "--output-on-failure", "-C", build_config),
                build_directory,
            )
        ]
    commands.extend(build_python_test_commands(root, python))
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
    platform_name: str | None = None,
) -> list[TestCommand]:
    """Build the selected executable when needed, then run CTest and Python tests."""
    if build_config not in {"Debug", "Release"}:
        raise ValueError(f"Unsupported CMake build configuration: {build_config!r}")
    root = Path(repository_root).resolve()
    executable = oneflow_executable_path(
        root, build_config, platform_name=platform_name
    )
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
