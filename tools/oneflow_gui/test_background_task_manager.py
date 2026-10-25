"""Tests for launching and stopping background task hosts."""

from __future__ import annotations

import sys
import tempfile
import time
import unittest
from pathlib import Path

from tools.oneflow_gui.background_task_manager import (
    request_background_stop,
    start_background_task,
)
from tools.oneflow_gui.background_task_registry import BackgroundTask, BackgroundTaskRegistry


class BackgroundTaskManagerTest(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.repository_root = Path(__file__).resolve().parents[2]
        self.registry = BackgroundTaskRegistry(self.root)
        self.log_file = self.root / "logs" / "task.log"
        self.registry.save(
            BackgroundTask.create(
                task_id="task-002",
                project_directory=self.root,
                executable=Path(sys.executable),
                execution_mode="source",
                log_file=self.log_file,
            )
        )

    def test_launch_detached_host_and_wait_for_terminal_state(self) -> None:
        start_background_task(
            root=self.root,
            task_id="task-002",
            command=[sys.executable, "-c", "print('managed')"],
            cwd=self.root,
            repository_root=self.repository_root,
        )
        deadline = time.monotonic() + 10
        while time.monotonic() < deadline:
            task = self.registry.load("task-002")
            if task.status in {"succeeded", "failed", "stopped", "lost"}:
                break
            time.sleep(0.05)
        self.assertEqual(self.registry.load("task-002").status, "succeeded")
        self.assertIn("managed", self.log_file.read_text(encoding="utf-8"))

    def test_stop_request_is_written_for_active_task(self) -> None:
        request_background_stop(self.root, "task-002")
        self.assertTrue((self.registry.records_directory / "task-002.stop").is_file())
        self.assertEqual(self.registry.load("task-002").status, "stopping")


if __name__ == "__main__":
    unittest.main()
