"""Tests for the serial case queue worker."""

from __future__ import annotations

import json
import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from tools.oneflow_gui.main_window import RunQueueWorker
from tools.oneflow_gui.models import CaseDefinition, CaseRunResult, RunState


class RunQueueWorkerTest(unittest.TestCase):
    def test_unexpected_runner_error_is_recorded_and_queue_continues(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            executable = root / "OneFLOW.exe"
            executable.write_bytes(b"placeholder")
            cases = [
                CaseDefinition("first", "First", root, root / "control.txt", None, None),
                CaseDefinition("second", "Second", root, root / "control.txt", None, None),
            ]
            worker = RunQueueWorker(cases, executable, root / "runs", "run-001")
            results = []
            states = []
            logs = []
            worker.case_result.connect(results.append)
            worker.case_state.connect(lambda case_id, state: states.append((case_id, state)))
            worker.log_line.connect(logs.append)

            def run_case_side_effect(case, *_args, **_kwargs):
                if case.case_id == "first":
                    raise RuntimeError("unexpected runner failure")
                return CaseRunResult(
                    case.case_id,
                    RunState.SUCCEEDED,
                    0,
                    root / "runs" / "run-001" / case.case_id,
                    "Process exited successfully.",
                )

            user_data = root / "user-data"
            with patch.dict(os.environ, {"ONEFLOW_HOME": str(user_data)}), patch(
                "tools.oneflow_gui.main_window.run_case", side_effect=run_case_side_effect
            ):
                worker.run()

            self.assertEqual(
                [(result.case_id, result.state) for result in results],
                [("first", RunState.FAILED), ("second", RunState.SUCCEEDED)],
            )
            self.assertIn(("first", RunState.FAILED.value), states)
            self.assertIn(("second", RunState.SUCCEEDED.value), states)
            self.assertTrue(any("unexpected runner failure" in line for line in logs))
            evidence_root = user_data / "workflow-data" / "executions" / "run-001"
            environment = json.loads(
                (evidence_root / "environment.json").read_text(encoding="utf-8")
            )
            manifest = json.loads(
                (evidence_root / "run-manifest.json").read_text(encoding="utf-8")
            )
            self.assertEqual(environment["capture_status"], "failed")
            self.assertEqual(manifest["status"], "failed")
            self.assertEqual(manifest["return_code"], 1)


if __name__ == "__main__":
    unittest.main()
