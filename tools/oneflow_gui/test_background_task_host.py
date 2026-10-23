"""Tests for detached command supervision."""

from __future__ import annotations

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from tools.oneflow_gui.background_task_host import run_host
from tools.oneflow_gui.background_task_registry import BackgroundTask, BackgroundTaskRegistry


class BackgroundTaskHostTest(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.registry = BackgroundTaskRegistry(self.root)
        self.log_file = self.root / "logs" / "task.log"
        self.registry.save(
            BackgroundTask.create(
                task_id="task-001",
                project_directory=self.root,
                executable=Path(sys.executable),
                execution_mode="source",
                log_file=self.log_file,
            )
        )

    def test_host_persists_success_and_captures_output(self) -> None:
        result = run_host(
            root=self.root,
            task_id="task-001",
            command=[sys.executable, "-c", "print('host-output')"],
            cwd=self.root,
            poll_interval=0.01,
        )

        task = self.registry.load("task-001")
        self.assertEqual(result, 0)
        self.assertEqual(task.status, "succeeded")
        self.assertIsNotNone(task.host_pid)
        self.assertIsNotNone(task.process_pid)
        self.assertIn("host-output", self.log_file.read_text(encoding="utf-8"))

    def test_host_records_nonzero_exit(self) -> None:
        result = run_host(
            root=self.root,
            task_id="task-001",
            command=[sys.executable, "-c", "raise SystemExit(7)"],
            cwd=self.root,
            poll_interval=0.01,
        )

        task = self.registry.load("task-001")
        self.assertEqual(result, 7)
        self.assertEqual(task.status, "failed")
        self.assertIn("code 7", task.message)


if __name__ == "__main__":
    unittest.main()
