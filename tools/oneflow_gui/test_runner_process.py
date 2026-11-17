"""Integration tests for case runner process-tree cancellation."""

from __future__ import annotations

import sys
import tempfile
import threading
import time
import unittest
from pathlib import Path
from unittest.mock import patch

from tools.oneflow_gui.models import CaseDefinition, RunState
from tools.oneflow_gui.runner import run_case


class CaseRunnerProcessTest(unittest.TestCase):
    def test_cancellation_stops_solver_descendant_processes(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "source-case"
            source.mkdir()
            (source / "control.txt").write_text("case input", encoding="utf-8")
            case = CaseDefinition(
                "process-tree-case",
                "Process Tree Case",
                source,
                source / "control.txt",
                None,
                None,
            )
            started_marker = root / "child-started.txt"
            survived_marker = root / "child-survived.txt"
            child_code = (
                "import pathlib, time; "
                f"pathlib.Path({str(started_marker)!r}).write_text('ready'); "
                "time.sleep(1); "
                f"pathlib.Path({str(survived_marker)!r}).write_text('alive')"
            )
            parent_code = (
                "import pathlib, subprocess, sys, time\n"
                f"subprocess.Popen([sys.executable, '-S', '-c', {child_code!r}])\n"
                f"started = pathlib.Path({str(started_marker)!r})\n"
                "deadline = time.monotonic() + 5\n"
                "while not started.exists() and time.monotonic() < deadline:\n"
                "    time.sleep(0.01)\n"
                "print('solver-started', flush=True)\n"
                "time.sleep(30)\n"
            )
            cancel_event = threading.Event()

            def cancel_after_child_starts() -> None:
                deadline = time.monotonic() + 5
                while not started_marker.is_file() and time.monotonic() < deadline:
                    time.sleep(0.01)
                if started_marker.is_file():
                    time.sleep(0.2)
                cancel_event.set()

            cancellation_thread = threading.Thread(
                target=cancel_after_child_starts,
                daemon=True,
            )
            cancellation_thread.start()
            with patch(
                "tools.oneflow_gui.runner.build_command",
                return_value=[sys.executable, "-S", "-c", parent_code],
            ):
                result = run_case(
                    case,
                    Path(sys.executable),
                    root / "workspace",
                    "run-001",
                    cancel_event,
                    lambda message: None,
                )
            cancellation_thread.join(timeout=5)

            self.assertEqual(result.state, RunState.CANCELLED)
            self.assertTrue(started_marker.is_file())
            time.sleep(1.2)
            self.assertFalse(survived_marker.exists())


if __name__ == "__main__":
    unittest.main()
