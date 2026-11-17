"""Tests for GUI-independent project execution."""

from __future__ import annotations

import json
import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from tools.oneflow_gui.background_project_runner import run_project
from tools.oneflow_gui.workflow_process import ProcessResult


class BackgroundProjectRunnerTest(unittest.TestCase):
    def test_runs_plain_project_without_qt(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            project = root / "case"
            project.mkdir()
            executable = root / "OneFLOW.exe"
            executable.touch()
            repository = root / "repository"
            repository.mkdir()

            user_data = root / "user-data"
            with patch.dict(os.environ, {"ONEFLOW_HOME": str(user_data)}), patch(
                "tools.oneflow_gui.background_project_runner.run_process",
                return_value=ProcessResult("succeeded", 0, "completed"),
            ) as run_process:
                result = run_project(
                    project,
                    executable,
                    "run-001",
                    "source",
                    None,
                    repository,
                )

            self.assertEqual(result, 0)
            run_process.assert_called_once()
            command = run_process.call_args.args[0]
            self.assertEqual(command, [str(executable.resolve()), "r", str(project.resolve())])
            self.assertFalse(run_process.call_args.kwargs["is_cancelled"]())
            evidence_root = user_data / "workflow-data" / "executions" / "run-001"
            environment = json.loads(
                (evidence_root / "environment.json").read_text(encoding="utf-8")
            )
            manifest = json.loads(
                (evidence_root / "run-manifest.json").read_text(encoding="utf-8")
            )
            self.assertEqual(environment["capture_status"], "failed")
            self.assertEqual(manifest["status"], "succeeded")
            self.assertEqual(manifest["return_code"], 0)


if __name__ == "__main__":
    unittest.main()
