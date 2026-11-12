"""Tests for cross-platform process-group lifecycle helpers."""

from __future__ import annotations

import os
import subprocess
import sys
import tempfile
import time
import unittest
from pathlib import Path
from unittest.mock import Mock

from tools.oneflow_gui.process_control import (
    process_group_options,
    terminate_and_reap,
    terminate_process_tree,
)


class ProcessControlTest(unittest.TestCase):
    def test_posix_process_group_options_start_a_new_session(self) -> None:
        self.assertEqual(
            process_group_options(platform_name="posix"),
            {"start_new_session": True},
        )

    def test_windows_process_group_options_create_a_new_process_group(self) -> None:
        self.assertEqual(
            process_group_options(platform_name="nt"),
            {"creationflags": 0x00000200},
        )

    def test_windows_termination_uses_taskkill_tree_mode(self) -> None:
        process = Mock()
        process.pid = 1234
        process.poll.return_value = None
        runner = Mock(return_value=subprocess.CompletedProcess([], 0))

        self.assertTrue(
            terminate_process_tree(
                process,
                platform_name="nt",
                command_runner=runner,
            )
        )
        self.assertEqual(
            runner.call_args.args[0],
            ["taskkill", "/PID", "1234", "/T", "/F"],
        )
        process.terminate.assert_not_called()

    def test_windows_termination_checks_tree_before_job_object(self) -> None:
        process = Mock()
        process.pid = 1234
        process.poll.return_value = None
        kernel32 = Mock()
        process._oneflow_job = (kernel32, 42)
        runner = Mock(return_value=subprocess.CompletedProcess([], 0))

        self.assertTrue(
            terminate_process_tree(
                process,
                platform_name="nt",
                command_runner=runner,
            )
        )

        self.assertEqual(
            runner.call_args.args[0],
            ["taskkill", "/PID", "1234", "/T", "/F"],
        )
        kernel32.TerminateJobObject.assert_not_called()

    def test_windows_termination_uses_job_object_when_taskkill_fails(self) -> None:
        process = Mock()
        process.pid = 1234
        process.poll.return_value = None
        kernel32 = Mock()
        kernel32.TerminateJobObject.return_value = True
        process._oneflow_job = (kernel32, 42)
        runner = Mock(return_value=subprocess.CompletedProcess([], 1))

        self.assertTrue(
            terminate_process_tree(
                process,
                platform_name="nt",
                command_runner=runner,
            )
        )

        kernel32.TerminateJobObject.assert_called_once_with(42, 1)
        process.terminate.assert_not_called()

    def test_windows_termination_falls_back_to_parent_process(self) -> None:
        process = Mock()
        process.pid = 1234
        process.poll.return_value = None
        runner = Mock(return_value=subprocess.CompletedProcess([], 1))

        self.assertTrue(
            terminate_process_tree(
                process,
                platform_name="windows",
                command_runner=runner,
            )
        )
        process.terminate.assert_called_once_with()

    def test_terminate_and_reap_waits_for_direct_child(self) -> None:
        process = Mock()
        process.poll.side_effect = [None, None, 0]
        process.wait.return_value = 0
        runner = Mock(return_value=subprocess.CompletedProcess([], 0))

        self.assertTrue(
            terminate_and_reap(
                process,
                platform_name="nt",
                command_runner=runner,
            )
        )

        runner.assert_called_once()
        process.wait.assert_called_once_with(timeout=3.0)
        process.kill.assert_not_called()

    def test_terminate_and_reap_escalates_after_timeout(self) -> None:
        process = Mock()
        process.poll.side_effect = [None, None, None, 0]
        process.wait.side_effect = [subprocess.TimeoutExpired("solver", 1.0), 0]
        runner = Mock(return_value=subprocess.CompletedProcess([], 0))

        self.assertTrue(
            terminate_and_reap(
                process,
                timeout=1.0,
                platform_name="nt",
                command_runner=runner,
            )
        )

        self.assertEqual(runner.call_count, 2)
        process.kill.assert_not_called()
        self.assertEqual(process.wait.call_count, 2)

    def test_terminate_and_reap_kills_child_if_tree_termination_does_not(self) -> None:
        process = Mock()
        process.poll.side_effect = [None, None, None, None, None, None]
        process.wait.side_effect = [subprocess.TimeoutExpired("solver", 1.0), 0]
        runner = Mock(return_value=subprocess.CompletedProcess([], 1))

        self.assertTrue(
            terminate_and_reap(
                process,
                timeout=1.0,
                platform_name="nt",
                command_runner=runner,
            )
        )

        process.kill.assert_called_once_with()
        process.wait.assert_called_with()

    def test_completed_process_is_not_terminated(self) -> None:
        process = Mock()
        process.poll.return_value = 0
        runner = Mock()

        self.assertFalse(
            terminate_process_tree(
                process,
                platform_name="nt",
                command_runner=runner,
            )
        )
        runner.assert_not_called()

    @unittest.skipUnless(os.name == "nt", "Requires Windows taskkill")
    def test_windows_termination_stops_descendant_processes(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            started_marker = Path(directory) / "child-started.txt"
            marker = Path(directory) / "child-survived.txt"
            release_marker = Path(directory) / "release-child.txt"
            child_code = (
                "import pathlib, time; "
                f"pathlib.Path({str(started_marker)!r}).write_text('ready'); "
                f"release = pathlib.Path({str(release_marker)!r})\n"
                "while not release.exists():\n"
                "    time.sleep(0.01)\n"
                f"pathlib.Path({str(marker)!r}).write_text('alive')"
            )
            parent_code = (
                "import pathlib, subprocess, sys, time\n"
                f"subprocess.Popen([sys.executable, '-S', '-c', {child_code!r}])\n"
                f"started = pathlib.Path({str(started_marker)!r})\n"
                "deadline = time.monotonic() + 5\n"
                "while not started.exists() and time.monotonic() < deadline:\n"
                "    time.sleep(0.01)\n"
                "print('child-started', flush=True)\n"
                "time.sleep(30)\n"
            )
            process = subprocess.Popen(
                [sys.executable, "-S", "-c", parent_code],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                **process_group_options(platform_name="nt"),
            )
            try:
                self.assertEqual(process.stdout.readline().strip(), "child-started")
                self.assertTrue(started_marker.is_file())
                self.assertTrue(
                    terminate_process_tree(process, platform_name="nt")
                )
                process.wait(timeout=10)
                # Release the child only after the tree-termination call. If
                # it survived, it must create the marker; no sleep-based race
                # can let the child finish before termination is attempted.
                release_marker.write_text("release", encoding="utf-8")
                deadline = time.monotonic() + 5
                while not marker.exists() and time.monotonic() < deadline:
                    time.sleep(0.01)
                self.assertFalse(marker.exists())
            finally:
                # Always release a surviving child, including when an earlier
                # assertion aborts the test before the normal release point.
                release_marker.write_text("release", encoding="utf-8")
                if process.poll() is None:
                    terminate_process_tree(process, platform_name="nt")
                    process.wait(timeout=10)
                if process.stdout is not None:
                    process.stdout.close()
                if process.stderr is not None:
                    process.stderr.close()

    @unittest.skipUnless(os.name == "posix", "Requires POSIX process groups")
    def test_posix_termination_stops_descendant_processes(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            started_marker = Path(directory) / "child-started.txt"
            marker = Path(directory) / "child-survived.txt"
            child_code = (
                "import pathlib, time; "
                f"pathlib.Path({str(started_marker)!r}).write_text('ready'); "
                "time.sleep(1); "
                f"pathlib.Path({str(marker)!r}).write_text('alive')"
            )
            parent_code = (
                "import pathlib, subprocess, sys, time\n"
                f"subprocess.Popen([sys.executable, '-c', {child_code!r}])\n"
                f"started = pathlib.Path({str(started_marker)!r})\n"
                "deadline = time.monotonic() + 5\n"
                "while not started.exists() and time.monotonic() < deadline:\n"
                "    time.sleep(0.01)\n"
                "print('child-started', flush=True)\n"
                "time.sleep(30)\n"
            )
            process = subprocess.Popen(
                [sys.executable, "-c", parent_code],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                start_new_session=True,
            )
            try:
                self.assertEqual(process.stdout.readline().strip(), "child-started")
                self.assertTrue(started_marker.is_file())
                self.assertTrue(
                    terminate_process_tree(process, platform_name="posix")
                )
                process.wait(timeout=5)
                time.sleep(1.2)
                self.assertFalse(marker.exists())
            finally:
                if process.poll() is None:
                    terminate_process_tree(process, platform_name="posix")
                    process.wait(timeout=5)
                if process.stdout is not None:
                    process.stdout.close()
                if process.stderr is not None:
                    process.stderr.close()


if __name__ == "__main__":
    unittest.main()
