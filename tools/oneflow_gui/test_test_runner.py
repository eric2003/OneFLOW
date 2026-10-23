"""Tests for building commands that reuse the repository's test runners."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from tools.oneflow_gui.test_runner import (
    available_suites,
    build_cmake_command,
    build_cmake_configure_command,
    build_default_test_commands,
    cmake_build_tree_is_configured,
    build_regression_command,
    build_unit_test_commands,
    stage_regression_tests,
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
        self.assertFalse(any("OneFLOW.exe" in " ".join(command.arguments) for command in commands))

    def test_configure_command_uses_build_directory_and_repository_root(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory).resolve()
            command = build_cmake_configure_command(root)
            self.assertEqual(command.working_directory, root / "build")
            self.assertEqual(
                command.arguments,
                (
                    "cmake", "..", "-D", "CMAKE_PREFIX_PATH=C:/dev/OneFLOW",
                    "-D", "CMAKE_BUILD_TYPE=Debug",
                ),
            )
            self.assertEqual(command.label, "CMake configure (Debug)")

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
            commands = build_default_test_commands(root, python="python-test")
            self.assertEqual(commands[0].label, "CMake build (Debug)")
            self.assertEqual(commands[1].label, "C++ unit tests (CTest)")

    def test_default_tests_do_not_rebuild_when_selected_executable_exists(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory).resolve()
            executable = root / "build" / "bin" / "Debug" / "OneFLOW.exe"
            executable.parent.mkdir(parents=True)
            executable.write_bytes(b"test")
            commands = build_default_test_commands(root, python="python-test")
            self.assertFalse(any(command.label.startswith("CMake ") for command in commands))
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

            staged_test_directory = root / "user-data" / "cfd-regression"
            staged_test_directory.mkdir(parents=True)
            (staged_test_directory / "suites").mkdir()
            (staged_test_directory / "suites" / "ci-serial.txt").write_text(
                "cfd/example\\n", encoding="utf-8"
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

    def test_regression_tests_are_copied_and_previous_workspace_is_replaced(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "repo"
            (root / "test" / "suites").mkdir(parents=True)
            (root / "test" / "suites" / "cpu-serial.txt").write_text(
                "cfd/example\\n", encoding="utf-8"
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


if __name__ == "__main__":
    unittest.main()
