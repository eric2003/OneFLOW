"""Tests for building commands that reuse the repository's test runners."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from tools.oneflow_gui.task_types import TaskType
from tools.oneflow_gui.test_runner import (
    available_suites,
    build_cmake_command,
    build_cmake_configure_command,
    build_ctest_test_commands,
    build_python_test_commands,
    build_default_test_commands,
    cmake_build_tree_is_configured,
    build_regression_command,
    build_unit_test_commands,
    format_command_arguments,
    format_final_test_status,
    oneflow_executable_path,
    stage_regression_tests,
    summarize_test_output,
)


class TestRunnerCommandTest(unittest.TestCase):
    def test_unit_tests_use_existing_ci_commands(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            commands = build_unit_test_commands(root, python="python-test")
        labels = [command.label for command in commands]
        self.assertIn("Mesh workflow tests", labels)
        self.assertIn("Workflow persistence tests", labels)
        self.assertIn("GUI test command builder tests", labels)
        self.assertIn("Runtime parameter consumer evidence tests", labels)
        python_modules = {
            command.arguments[-1]
            for command in commands
            if len(command.arguments) >= 4 and command.arguments[1:3] == ("-m", "unittest")
        }
        self.assertTrue({
            "tools.verification.test_curve_catalog",
            "tools.oneflow_gui.test_parameter_consumer_catalog",
            "tools.verification.test_contract_catalog",
            "tools.oneflow_gui.test_mesh_validation",
            "grid.tests.python.stl_to_gmsh.test_stl_to_gmsh_geometry",
            "tools.verification.test_cmake_test_catalog",
            "tools.oneflow_gui.test_contract_workbench_dialog",
            "tools.oneflow_gui.test_source_editor",
            "tools.oneflow_gui.test_test_inventory",
            "tools.oneflow_gui.test_test_diagnostics",
            "tools.oneflow_gui.test_executable_discovery",
            "tools.oneflow_gui.test_process_control",
            "tools.oneflow_gui.test_runner_process",
        }.issubset(python_modules))
        self.assertTrue(all(command.task_type is TaskType.UNIT_TESTS for command in commands))
        self.assertFalse(any("OneFLOW.exe" in " ".join(command.arguments) for command in commands))

    def test_ctest_and_python_commands_can_be_built_independently(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory).resolve()
            build = root / "build"
            build.mkdir()
            (build / "CTestTestfile.cmake").write_text("", encoding="utf-8")

            ctest_commands = build_ctest_test_commands(root, build_config="Release")
            python_commands = build_python_test_commands(root, python="python-test")

            self.assertEqual(len(ctest_commands), 1)
            self.assertEqual(ctest_commands[0].arguments, (
                "ctest", "--output-on-failure", "-C", "Release"
            ))
            self.assertEqual(ctest_commands[0].working_directory, build)
            self.assertTrue(python_commands)
            self.assertTrue(all(command.arguments[0] == "python-test" for command in python_commands))
            self.assertFalse(any(command.arguments[0] == "ctest" for command in python_commands))

    def test_ctest_builder_reports_unconfigured_build_tree(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory).resolve()
            self.assertEqual(build_ctest_test_commands(root), [])
            self.assertTrue(build_python_test_commands(root))

    def test_configure_command_uses_build_directory_and_repository_root(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory).resolve()
            command = build_cmake_configure_command(root)
            self.assertEqual(command.working_directory, root / "build")
            self.assertEqual(
                command.arguments,
                (
                    "cmake", "..", "-D", "CMAKE_BUILD_TYPE=Debug",
                ),
            )
            self.assertEqual(command.label, "CMake configure (Debug)")

    def test_linux_executable_path_uses_single_configuration_output(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory).resolve()
            self.assertEqual(
                oneflow_executable_path(root, "Debug", platform_name="posix"),
                root / "build" / "bin" / "OneFLOW",
            )
            executable = oneflow_executable_path(root, "Release", platform_name="posix")
            executable.parent.mkdir(parents=True)
            executable.write_bytes(b"test")
            commands = build_default_test_commands(
                root, python="python-test", build_config="Release", platform_name="posix"
            )
            self.assertFalse(any(command.label.startswith("CMake build") for command in commands))

    def test_posix_multi_configuration_generator_uses_config_subdirectory(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory).resolve()
            build = root / "build"
            build.mkdir()
            (build / "CMakeCache.txt").write_text(
                "CMAKE_CONFIGURATION_TYPES:STRING=Debug;Release\n",
                encoding="utf-8",
            )
            expected = build / "bin" / "Release" / "OneFLOW"
            self.assertEqual(
                oneflow_executable_path(root, "Release", platform_name="posix"),
                expected,
            )

    def test_command_log_uses_platform_appropriate_quoting(self) -> None:
        arguments = ("python", "path with spaces/test.py", "a'b")
        self.assertEqual(
            format_command_arguments(arguments, platform_name="posix"),
            "python 'path with spaces/test.py' 'a'\"'\"'b'",
        )
        self.assertEqual(
            format_command_arguments(("tool.exe", "path with spaces"), platform_name="nt"),
            'tool.exe "path with spaces"',
        )

    def test_release_configuration_is_used_for_cmake_configure(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            command = build_cmake_configure_command(Path(directory), "Release")
            self.assertEqual(command.label, "CMake configure (Release)")
            self.assertEqual(command.arguments[-1], "CMAKE_BUILD_TYPE=Release")

    def test_build_tree_requires_cache_and_generated_build_file(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            build = Path(directory) / "build"
            build.mkdir()
            self.assertFalse(cmake_build_tree_is_configured(build))
            (build / "CMakeCache.txt").write_text("", encoding="utf-8")
            self.assertFalse(cmake_build_tree_is_configured(build))
            (build / "build.ninja").write_text("", encoding="utf-8")
            self.assertFalse(cmake_build_tree_is_configured(build))
            (build / "CTestTestfile.cmake").write_text("", encoding="utf-8")
            self.assertTrue(cmake_build_tree_is_configured(build))

    def test_default_tests_configure_build_and_run_ctest_when_executable_is_missing(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory).resolve()
            commands = build_default_test_commands(root, python="python-test")
            labels = [command.label for command in commands]
            self.assertEqual(labels[0:3], [
                "CMake configure (Debug)",
                "CMake build (Debug)",
                "C++ unit tests (CTest)",
            ])
            self.assertIn("GUI test command builder tests", labels)
            self.assertEqual(commands[0].working_directory, root / "build")
            self.assertEqual(commands[1].working_directory, root)
            self.assertEqual(commands[2].working_directory, root / "build")
            self.assertEqual(commands[2].arguments, (
                "ctest", "--output-on-failure", "-C", "Debug"
            ))

    def test_default_tests_skip_configure_but_build_when_build_tree_exists(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory).resolve()
            build = root / "build"
            build.mkdir()
            (build / "CMakeCache.txt").write_text("", encoding="utf-8")
            (build / "build.ninja").write_text("", encoding="utf-8")
            (build / "CTestTestfile.cmake").write_text("", encoding="utf-8")
            commands = build_default_test_commands(root, python="python-test")
            self.assertEqual(commands[0].label, "CMake build (Debug)")
            self.assertEqual(commands[1].label, "C++ unit tests (CTest)")

    def test_default_tests_do_not_rebuild_when_selected_executable_exists(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory).resolve()
            executable = oneflow_executable_path(root, "Debug", platform_name="nt")
            executable.parent.mkdir(parents=True)
            executable.write_bytes(b"test")
            commands = build_default_test_commands(root, python="python-test", platform_name="nt")
            build_action_prefixes = (
                "CMake configure (",
                "CMake build (",
                "CMake rebuild (",
            )
            self.assertFalse(
                any(command.label.startswith(build_action_prefixes) for command in commands)
            )
            self.assertFalse(any(command.label.startswith("C++") for command in commands))

    def test_ctest_is_included_only_when_build_is_configured(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            build = root / "build"
            build.mkdir()
            self.assertFalse(any(command.label.startswith("C++") for command in build_unit_test_commands(root)))
            (build / "CTestTestfile.cmake").write_text("", encoding="utf-8")
            commands = build_unit_test_commands(root)
            ctest = next(command for command in commands if command.label.startswith("C++"))
            self.assertEqual(ctest.working_directory, build.resolve())
            self.assertEqual(
                ctest.arguments,
                ("ctest", "--output-on-failure", "-C", "Debug"),
            )

            release_commands = build_unit_test_commands(
                root, build_config="Release"
            )
            release_ctest = next(
                command for command in release_commands
                if command.label.startswith("C++")
            )
            self.assertEqual(
                release_ctest.arguments,
                ("ctest", "--output-on-failure", "-C", "Release"),
            )
            self.assertEqual(release_ctest.working_directory, build.resolve())

    def test_ctest_rejects_unknown_build_configuration(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaisesRegex(ValueError, "Unsupported CMake build configuration"):
                build_unit_test_commands(Path(directory), build_config="RelWithDebInfo")

    def test_cmake_build_and_rebuild_use_selected_configuration(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            build = build_cmake_command(root, "Debug")
            self.assertEqual(
                build.arguments,
                ("cmake", "--build", str((root / "build").resolve()), "--config", "Debug"),
            )
            self.assertEqual(build.working_directory, root.resolve())
            self.assertEqual(build.label, "CMake build (Debug)")

            rebuild = build_cmake_command(root, "Release", rebuild=True)
            self.assertEqual(
                rebuild.arguments,
                (
                    "cmake", "--build", str((root / "build").resolve()),
                    "--config", "Release", "--clean-first",
                ),
            )
            self.assertEqual(rebuild.label, "CMake rebuild (Release)")

    def test_cmake_build_rejects_unknown_configuration(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaisesRegex(ValueError, "Unsupported CMake build configuration"):
                build_cmake_command(Path(directory), "RelWithDebInfo")

    def test_available_suites_comes_from_existing_manifests(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            suites = Path(directory) / "test" / "suites"
            suites.mkdir(parents=True)
            (suites / "ci-serial.txt").write_text("@cpu-serial.txt\n", encoding="utf-8")
            (suites / "cpu-serial.txt").write_text("cfd/example\n", encoding="utf-8")
            self.assertEqual(available_suites(Path(directory)), ["ci-serial.txt", "cpu-serial.txt"])

    def test_regression_command_uses_test_py_and_existing_suite(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            test_dir = root / "test"
            (test_dir / "suites").mkdir(parents=True)
            (test_dir / "suites" / "ci-serial.txt").write_text("cfd/example\n", encoding="utf-8")
            executable = root / "Program Files" / "OneFLOW.exe"
            executable.parent.mkdir(parents=True)
            executable.write_bytes(b"test")
            command = build_regression_command(
                root, executable, "ci-serial.txt", python="python-test", platform_name="nt"
            )
            self.assertEqual(command.working_directory, test_dir.resolve())
            self.assertEqual(command.arguments[0:2], ("python-test", "test.py"))
            self.assertEqual(command.arguments[2], "mpiexec -n 1")
            self.assertEqual(command.arguments[3], '"' + str(executable.resolve()) + '"')
            self.assertEqual(command.arguments[4], "suites/ci-serial.txt")
            self.assertIs(command.task_type, TaskType.CFD_GRID_REGRESSION)

            staged_test_directory = root / "user-data" / "cfd-regression"
            staged_test_directory.mkdir(parents=True)
            (staged_test_directory / "suites").mkdir()
            (staged_test_directory / "suites" / "ci-serial.txt").write_text(
                "cfd/example\n", encoding="utf-8"
            )
            staged_command = build_regression_command(
                root,
                executable,
                "ci-serial.txt",
                python="python-test",
                platform_name="nt",
                test_directory=staged_test_directory,
            )
            self.assertEqual(staged_command.working_directory, staged_test_directory.resolve())
            self.assertEqual(staged_command.arguments[4], "suites/ci-serial.txt")

    def test_regression_command_uses_posix_mpi_and_shell_quoting(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            test_dir = root / "test"
            (test_dir / "suites").mkdir(parents=True)
            (test_dir / "suites" / "linux.txt").write_text("cfd/example\n", encoding="utf-8")
            executable = root / "solver with spaces" / "OneFLOW"
            executable.parent.mkdir(parents=True)
            executable.write_bytes(b"test")

            command = build_regression_command(
                root,
                executable,
                "linux.txt",
                python="python3",
                platform_name="posix",
            )

            self.assertEqual(command.arguments[2], "mpirun -np 1")
            self.assertEqual(command.arguments[3], "'" + str(executable.resolve()) + "'")
            self.assertEqual(command.arguments[4], "suites/linux.txt")

    def test_regression_tests_are_copied_and_previous_workspace_is_replaced(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "repo"
            (root / "test" / "suites").mkdir(parents=True)
            (root / "test" / "suites" / "cpu-serial.txt").write_text(
                "cfd/example\n", encoding="utf-8"
            )
            (root / "test" / "cfd" / "example").mkdir(parents=True)
            (root / "test" / "cfd" / "example" / "config.json").write_text(
                "{}", encoding="utf-8"
            )
            workspace = Path(directory) / "user-data" / "cfd-regression"
            workspace.mkdir(parents=True)
            (workspace / "stale-output.dat").write_text("stale", encoding="utf-8")

            staged = stage_regression_tests(root, workspace)

            self.assertEqual(staged, workspace.resolve())
            self.assertFalse((staged / "stale-output.dat").exists())
            self.assertTrue((staged / "suites" / "cpu-serial.txt").is_file())
            self.assertTrue((staged / "cfd" / "example" / "config.json").is_file())
            self.assertFalse((root / "test" / "stale-output.dat").exists())

    def test_regression_workspace_must_be_outside_source_repository(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "repo"
            (root / "test").mkdir(parents=True)
            with self.assertRaisesRegex(ValueError, "outside the source repository"):
                stage_regression_tests(root, root / "test" / "run-copy")

    def test_regression_rejects_unknown_suite_and_missing_executable(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "test" / "suites").mkdir(parents=True)
            (root / "test" / "suites" / "ci-serial.txt").write_text(
                "cfd/example\n", encoding="utf-8"
            )
            with self.assertRaisesRegex(ValueError, "Unknown regression suite"):
                build_regression_command(root, root / "missing.exe", "missing.txt")
            executable = root / "OneFLOW.exe"
            executable.write_bytes(b"test")
            with self.assertRaisesRegex(ValueError, "Unknown regression suite"):
                build_regression_command(root, executable, "../outside.txt")
            deferred = build_regression_command(
                root,
                root / "build" / "bin" / "Debug" / "OneFLOW.exe",
                "ci-serial.txt",
                validate_executable=False,
            )
            self.assertEqual(deferred.arguments[1], "test.py")


    def test_unittest_output_reports_case_counts_and_failure_traceback(self) -> None:
        report = summarize_test_output([
            "F..",
            "FAIL: test_case_domain_displays_typed_case_setup_tab (sample.TestCase)",
            "Traceback (most recent call last):",
            '  File "test_sample.py", line 10, in test_case',
            "AssertionError: expected label was missing",
            "Ran 3 tests in 0.153s",
            "FAILED (failures=1)",
        ])
        self.assertIn("Python unittest cases: 3 run, 2 passed, 1 failed, 0 skipped.", report)
        self.assertIn("Failed unittest cases:", report)
        self.assertIn(
            "  - test_case_domain_displays_typed_case_setup_tab (sample.TestCase)",
            report,
        )
        self.assertIn("Failure details:", report)
        self.assertIn("AssertionError: expected label was missing", report)

    def test_unittest_output_lists_every_failed_case_name(self) -> None:
        report = summarize_test_output([
            "FAIL: test_first (sample.FirstTest)",
            "ERROR: test_second (sample.SecondTest)",
            "Ran 5 tests in 0.1s",
            "FAILED (failures=1, errors=1)",
        ])
        self.assertIn("2 failed", report[0])
        self.assertIn("  - test_first (sample.FirstTest)", report)
        self.assertIn("  - test_second (sample.SecondTest)", report)

    def test_ctest_output_reports_passed_and_failed_test_counts(self) -> None:
        report = summarize_test_output([
            "67% tests passed, 1 tests failed out of 3",
            "The following tests FAILED:",
            "  2 - FaceTopologyTest (Failed)",
        ])
        self.assertIn("CTest cases: 3 run, 2 passed, 1 failed, 0 skipped.", report)
        self.assertIn("Failure details:", report)
        self.assertIn("  2 - FaceTopologyTest (Failed)", report)

    def test_ctest_output_reports_skipped_test_counts(self) -> None:
        report = summarize_test_output([
            "100% tests passed, 0 tests failed out of 3",
            "The following tests did not run:",
            "  3 - OptionalHardwareTest (Skipped)",
        ])
        self.assertIn("CTest cases: 3 run, 2 passed, 0 failed, 1 skipped.", report)

    def test_final_status_includes_command_and_framework_case_counts(self) -> None:
        report = "\n".join([
            "Test command summary:",
            "  Total: 8",
            "  Passed: 5",
            "  Failed: 1",
            "  Not run: 2",
            "Python unittest cases: 26 run, 24 passed, 1 failed, 1 skipped.",
            "CTest cases: 10 run, 9 passed, 0 failed, 1 skipped.",
        ])
        status = format_final_test_status(False, report)
        self.assertEqual(
            status,
            "FINAL STATUS: FAIL | Commands: 5 passed, 1 failed, 2 not run (8 total) "
            "| Python unittest: 24 passed, 1 failed, 1 skipped (26 total) "
            "| CTest: 9 passed, 0 failed, 1 skipped (10 total)",
        )

    def test_final_status_without_case_summary_keeps_command_counts(self) -> None:
        status = format_final_test_status(
            True,
            "Test command summary:\n  Total: 2\n  Passed: 2\n"
            "  Failed: 0\n  Not run: 0",
        )
        self.assertEqual(
            status,
            "FINAL STATUS: PASS | Commands: 2 passed, 0 failed, 0 not run (2 total)",
        )

    def test_multiple_unittest_runs_aggregate_each_failure_count(self) -> None:
        report = summarize_test_output([
            "Ran 10 tests in 0.2s",
            "OK",
            "Ran 21 tests in 0.1s",
            "FAILED (failures=2)",
        ])
        self.assertIn(
            "Python unittest cases: 31 run, 29 passed, 2 failed, 0 skipped.",
            report,
        )

    def test_successful_unittest_output_reports_case_count_without_failure(self) -> None:
        report = summarize_test_output([
            "Ran 12 tests in 0.8s",
            "OK",
        ])
        self.assertEqual(report, ("Python unittest cases: 12 run, 12 passed, 0 failed, 0 skipped.",))


if __name__ == "__main__":
    unittest.main()
