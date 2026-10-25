"""Launch and control detached background task hosts."""

from __future__ import annotations

import json
import os
import subprocess
import sys
import time
from dataclasses import replace
from pathlib import Path
from typing import Sequence

from .background_task_registry import ACTIVE_STATUSES, BackgroundTaskRegistry, utc_now


def start_background_task(
    *,
    root: Path,
    task_id: str,
    command: Sequence[str],
    cwd: Path,
    repository_root: Path,
    python_executable: str | None = None,
    startup_timeout: float = 5.0,
) -> int:
    """Launch a detached host and wait until its registry record acknowledges it."""
    if not command or not all(isinstance(item, str) for item in command):
        raise ValueError("Background command must be a non-empty string sequence.")
    registry = BackgroundTaskRegistry(root)
    task = registry.load(task_id)
    if task.status not in ACTIVE_STATUSES:
        raise ValueError(f"Task {task_id!r} is not in an active state.")
    config_directory = root / "background-task-config"
    config_directory.mkdir(parents=True, exist_ok=True)
    config_file = config_directory / f"{task_id}.json"
    config_file.write_text(
        json.dumps({"command": list(command)}, ensure_ascii=True, indent=2) + "\n",
        encoding="utf-8",
    )
    python = python_executable or sys.executable
    host_command = [
        python,
        "-m",
        "tools.oneflow_gui.background_task_host",
        "--root",
        str(root.resolve()),
        "--task-id",
        task_id,
        "--command-file",
        str(config_file.resolve()),
        "--cwd",
        str(cwd.resolve()),
    ]
    options = {
        "stdin": subprocess.DEVNULL,
        "stdout": subprocess.DEVNULL,
        "stderr": subprocess.DEVNULL,
        "cwd": str(repository_root.resolve()),
        "close_fds": True,
    }
    if os.name == "nt":
        options["creationflags"] = (
            subprocess.DETACHED_PROCESS | subprocess.CREATE_NEW_PROCESS_GROUP
        )
    else:
        options["start_new_session"] = True
    process = subprocess.Popen(host_command, **options)
    deadline = time.monotonic() + startup_timeout
    while time.monotonic() < deadline:
        current = registry.load(task_id)
        if current.host_pid == process.pid or current.status not in ACTIVE_STATUSES:
            if current.host_pid == process.pid or current.status in {
                "succeeded",
                "failed",
                "stopped",
            }:
                return process.pid
            break
        time.sleep(0.05)
    try:
        process.poll()
    except OSError:
        pass
    current = registry.load(task_id)
    if current.host_pid != process.pid and current.status in ACTIVE_STATUSES:
        raise RuntimeError(
            f"Background host did not acknowledge task {task_id!r} "
            f"within {startup_timeout:.1f} seconds."
        )
    return process.pid


def request_background_stop(root: Path, task_id: str) -> None:
    """Ask the host to stop its active process using a persistent request file."""
    registry = BackgroundTaskRegistry(root)
    task = registry.load(task_id)
    if task.status not in ACTIVE_STATUSES:
        return
    stop_file = registry.records_directory / f"{task_id}.stop"
    stop_file.write_text(utc_now() + "\n", encoding="utf-8")
    registry.save(
        replace(
            task,
            status="stopping",
            updated_at=utc_now(),
            message="Stop requested by the user.",
        )
    )
