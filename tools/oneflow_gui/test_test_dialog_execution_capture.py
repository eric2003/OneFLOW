"""Integration tests for execution evidence at the test worker boundary."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

from tools.oneflow_gui.test_dialog import TestCommandWorker
from tools.oneflow_gui.test_runner import TestCommand


class TestCommandWorkerEvidenceTest(unittest.TestCase):
    def _command(self, root: Path) -> TestCommand:
        return TestCommand(
            "Sample unittest suite",
            ("python", "-m", "unittest", "tools.sample_tests"),
            root,
        )

    def _process(self, return_code: int = 0) -> MagicMock:
        process = MagicMock()
        process.stdout = iter(("Ran 1 test in 0.01s\n", "OK\n"))
        process.wait.return_value = return_code
        process.__enter__.return_value = process
        process.__exit__.return_value = False
        return process

    def test_worker_persists_evidence_after_unittest_process_completes(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            process = self._process()
            paths = SimpleNamespace(
                verification_execution_records=root / ".oneflow" / "verification" / "execution-records.jsonl",
                test_run_logs=root / ".oneflow" / "test-runs" / "logs",
            )
            with (
                patch(
                    "tools.oneflow_gui.test_dialog.subprocess.Popen",
                    return_value=process,
                ),
                patch(
                    "tools.oneflow_gui.test_dialog.UserDataPaths.discover",
                    return_value=paths,
                ),
                patch(
                    "tools.oneflow_gui.test_dialog.capture_python_unittest_execution",
                    return_value=SimpleNamespace(
                        record_id="record-1", outcome="passed"
                    ),
                ) as capture,
            ):
                worker = TestCommandWorker([self._command(root)], build_config="Debug")
                worker.run()

            capture.assert_called_once()
            arguments = capture.call_args.kwargs
            self.assertEqual(arguments["exit_code"], 0)
            self.assertEqual(arguments["output_lines"], [
                "Ran 1 test in 0.01s",
                "OK",
            ])
            self.assertEqual(arguments["configuration"], "Debug")
            self.assertEqual(
                arguments["artifact_directory"], paths.test_run_logs
            )
            self.assertFalse(arguments["cancelled"])

    def test_process_launch_failure_does_not_create_execution_evidence(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            with (
                patch(
                    "tools.oneflow_gui.test_dialog.subprocess.Popen",
                    side_effect=OSError("python executable unavailable"),
                ),
                patch(
                    "tools.oneflow_gui.test_dialog.capture_python_unittest_execution"
                ) as capture,
            ):
                worker = TestCommandWorker([self._command(root)])
                worker.run()

            capture.assert_not_called()

    def test_cancellation_before_launch_does_not_create_execution_evidence(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            with (
                patch("tools.oneflow_gui.test_dialog.subprocess.Popen") as popen,
                patch(
                    "tools.oneflow_gui.test_dialog.capture_python_unittest_execution"
                ) as capture,
            ):
                worker = TestCommandWorker([self._command(root)])
                worker.request_stop()
                worker.run()

            popen.assert_not_called()
            capture.assert_not_called()


if __name__ == "__main__":
    unittest.main()
