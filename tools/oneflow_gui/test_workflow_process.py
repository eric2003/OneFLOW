"""Tests for Qt-independent external process execution."""

from __future__ import annotations

import sys
import tempfile
import time
import unittest
from pathlib import Path
from unittest.mock import MagicMock, patch

from tools.oneflow_gui.workflow_process import run_process


class WorkflowProcessTest(unittest.TestCase):
    def setUp(self) -> None:
        self._temporary_directory = tempfile.TemporaryDirectory()
        self.addCleanup(self._temporary_directory.cleanup)
        self._cwd = Path(self._temporary_directory.name)

    def _run(self, code: str, cancelled=lambda: False):
        output: list[str] = []
        result = run_process(
            [sys.executable, "-c", code],
            cwd=self._cwd,
            on_output=output.append,
            is_cancelled=cancelled,
        )
        return result, output

    def test_success_streams_output_and_reports_success(self) -> None:
        result, output = self._run("print('mesh ready')")
        self.assertEqual(result.status, "succeeded")
        self.assertEqual(result.return_code, 0)
        self.assertIn("mesh ready", output)

    def test_nonzero_exit_is_a_failed_result(self) -> None:
        result, output = self._run("print('solver failed'); raise SystemExit(7)")
        self.assertEqual(result.status, "failed")
        self.assertEqual(result.return_code, 7)
        self.assertIn("solver failed", output)

    def test_cancellation_terminates_process(self) -> None:
        start = time.monotonic()
        result, _ = self._run(
            "import time; print('started', flush=True); time.sleep(10)",
            cancelled=lambda: time.monotonic() - start > 0.2,
        )
        self.assertEqual(result.status, "cancelled")
        self.assertNotEqual(result.return_code, 0)

    def test_windows_process_suppresses_console_window(self) -> None:
        process = MagicMock()
        process.stdout = iter([])
        process.poll.return_value = 0
        process.returncode = 0
        with (
            patch("tools.oneflow_gui.workflow_process.os.name", "nt"),
            patch(
                "tools.oneflow_gui.workflow_process.subprocess.CREATE_NO_WINDOW",
                0x08000000,
                create=True,
            ),
            patch(
                "tools.oneflow_gui.workflow_process.subprocess.Popen",
                return_value=process,
            ) as popen,
        ):
            result = run_process(
                ["OneFLOW.exe", "r", "case"],
                cwd=self._cwd,
                on_output=lambda line: None,
                is_cancelled=lambda: False,
            )

        self.assertEqual(result.status, "succeeded")
        self.assertEqual(popen.call_args.kwargs["creationflags"], 0x08000000)

    def test_empty_command_is_rejected(self) -> None:
        with self.assertRaisesRegex(ValueError, "non-empty string arguments"):
            run_process(
                [],
                cwd=self._cwd,
                on_output=lambda line: None,
                is_cancelled=lambda: False,
            )


if __name__ == "__main__":
    unittest.main()
