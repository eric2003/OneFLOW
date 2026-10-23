"""Tests for GUI-independent project execution."""

from __future__ import annotations

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

            with patch(
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


if __name__ == "__main__":
    unittest.main()
