"""Tests for real Python unittest execution evidence capture."""

from __future__ import annotations

import hashlib
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from tools.oneflow_gui.execution_capture import (
    capture_python_unittest_execution,
    is_python_unittest_command,
)
from tools.oneflow_gui.test_runner import TestCommand


class PythonExecutionCaptureTest(unittest.TestCase):
    def _command(self, root: Path) -> TestCommand:
        return TestCommand(
            "Verification execution record tests",
            (
                "python",
                "-m",
                "unittest",
                "tools.verification.test_execution_records",
            ),
            root,
        )

    def test_captures_completed_unittest_command_and_log_artifact(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            command = self._command(root)
            records_path = root / ".oneflow" / "records.jsonl"
            artifacts = root / ".oneflow" / "logs"
            git_values = {
                ("remote", "get-url", "origin"): "https://github.com/eric2003/OneFLOW.git",
                ("rev-parse", "HEAD"): "abc123",
                ("branch", "--show-current"): "hexin",
            }

            def git_result(arguments, **_kwargs):
                value = git_values.get(tuple(arguments[1:]), "")
                result = type("Result", (), {})()
                result.stdout = value
                result.returncode = 0 if value else 1
                return result

            with patch(
                "tools.oneflow_gui.execution_capture.subprocess.run",
                side_effect=git_result,
            ):
                record = capture_python_unittest_execution(
                    command,
                    repository_root=root,
                    records_path=records_path,
                    artifact_directory=artifacts,
                    started_at="2026-11-02T01:00:00+08:00",
                    finished_at="2026-11-02T01:00:02+08:00",
                    exit_code=0,
                    output_lines=("Ran 12 tests in 0.1s", "OK"),
                    configuration="Debug",
                )

            self.assertEqual(record.outcome, "passed")
            self.assertEqual(record.framework, "unittest")
            self.assertEqual(record.source_path, "tools/verification/test_execution_records.py")
            self.assertEqual(record.repository, "eric2003/OneFLOW")
            self.assertEqual(record.commit, "abc123")
            self.assertEqual(record.branch, "hexin")
            self.assertEqual(record.configuration, "Debug")
            self.assertEqual(record.contract_ids, ())
            self.assertEqual(record.command, command.arguments)
            artifact = record.artifacts[0]
            log_path = root / artifact.path
            self.assertTrue(log_path.is_file())
            self.assertEqual(
                artifact.sha256,
                hashlib.sha256(log_path.read_bytes()).hexdigest(),
            )
            self.assertIn("Ran 12 tests", log_path.read_text(encoding="utf-8"))

    def test_records_failed_and_cancelled_process_outcomes(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            command = self._command(root)
            common = {
                "repository_root": root,
                "records_path": root / "records.jsonl",
                "artifact_directory": root / "logs",
                "started_at": "2026-11-02T01:00:00+08:00",
                "finished_at": "2026-11-02T01:00:01+08:00",
                "output_lines": ("FAILED",),
            }

            failed = capture_python_unittest_execution(
                command, exit_code=1, **common
            )
            cancelled = capture_python_unittest_execution(
                command, exit_code=-15, cancelled=True, **common
            )

            self.assertEqual(failed.outcome, "failed")
            self.assertEqual(cancelled.outcome, "cancelled")

    def test_only_explicit_unittest_commands_are_captured(self) -> None:
        unittest_command = self._command(Path("."))
        ctest_command = TestCommand(
            "C++ unit tests (CTest)",
            ("ctest", "--output-on-failure"),
            Path("."),
        )
        script_command = TestCommand(
            "Python script",
            ("python", "tools/run_helper.py"),
            Path("."),
        )

        self.assertTrue(is_python_unittest_command(unittest_command))
        self.assertFalse(is_python_unittest_command(ctest_command))
        self.assertFalse(is_python_unittest_command(script_command))


if __name__ == "__main__":
    unittest.main()
