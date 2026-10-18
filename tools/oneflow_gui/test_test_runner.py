"""Tests for building commands that reuse the repository's test runners."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from tools.oneflow_gui.test_runner import (
    available_suites,
    build_regression_command,
    build_unit_test_commands,
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
            self.assertEqual(ctest.arguments[:2], ("ctest", "--output-on-failure"))

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
            executable = root / "OneFLOW.exe"
            executable.write_bytes(b"test")
            command = build_regression_command(
                root, executable, "ci-serial.txt", python="python-test", platform_name="nt"
            )
            self.assertEqual(command.working_directory, test_dir.resolve())
            self.assertEqual(command.arguments[0:2], ("python-test", "test.py"))
            self.assertEqual(command.arguments[2], "mpiexec -n 1")
            self.assertEqual(command.arguments[4], "suites/ci-serial.txt")

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


if __name__ == "__main__":
    unittest.main()
