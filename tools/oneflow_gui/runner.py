"""Process execution service for isolated OneFLOW case workspaces."""

from __future__ import annotations

import queue
import shutil
import subprocess
import threading
from datetime import datetime
from pathlib import Path
from typing import Callable
from uuid import uuid4

from .models import CaseDefinition, CaseRunResult, RunState


LogCallback = Callable[[str], None]


def build_command(executable: Path, case_directory: Path) -> list[str]:
    """Build the standard OneFLOW run command without shell interpolation."""

    return [str(executable.resolve()), "r", str(case_directory.resolve())]


def create_run_id() -> str:
    """Return a sortable, collision-resistant identifier for one batch."""

    timestamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    return f"{timestamp}-{uuid4().hex[:8]}"


def _copy_case_to_workspace(
    case: CaseDefinition,
    workspace_root: Path,
    run_id: str,
) -> Path:
    """Copy source inputs so generated files never overwrite checked-in cases."""

    destination = workspace_root / run_id / case.case_id
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists():
        raise FileExistsError(f"Case workspace already exists: {destination}")

    shutil.copytree(
        case.case_directory,
        destination,
        ignore=shutil.ignore_patterns("__pycache__", ".git"),
    )
    return destination


def run_case(
    case: CaseDefinition,
    executable: Path,
    workspace_root: Path,
    run_id: str,
    cancel_event: threading.Event,
    log: LogCallback,
    *,
    execution_mode: str = "default_workspace",
) -> CaseRunResult:
    """Run one case, forward combined output, and support cooperative stopping."""

    executable = executable.expanduser().resolve()
    if not executable.is_file():
        return CaseRunResult(
            case.case_id, RunState.FAILED, None, workspace_root,
            f"OneFLOW executable not found: {executable}",
        )
    if cancel_event.is_set():
        return CaseRunResult(
            case.case_id, RunState.CANCELLED, None, workspace_root,
            "Cancelled before process launch.",
        )

    try:
        if execution_mode == "source":
            case_workspace = case.case_directory.resolve()
            if not case_workspace.is_dir():
                raise NotADirectoryError(f"Source case directory does not exist: {case_workspace}")
        elif execution_mode in {"default_workspace", "custom_workspace"}:
            case_workspace = _copy_case_to_workspace(case, workspace_root, run_id)
        else:
            raise ValueError(f"Unsupported case execution mode: {execution_mode!r}")
        command = build_command(executable, case_workspace)
        log(f"[{case.case_id}] Command: {' '.join(command)}")
        log(f"[{case.case_id}] Working directory: {executable.parent}")

        process = subprocess.Popen(
            command,
            cwd=executable.parent,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            errors="replace",
            bufsize=1,
        )
        output_queue: queue.Queue[str | None] = queue.Queue()

        def read_output() -> None:
            """Drain the pipe on a dedicated thread to avoid blocking cancellation."""
            assert process.stdout is not None
            try:
                for line in process.stdout:
                    output_queue.put(line.rstrip("\r\n"))
            finally:
                output_queue.put(None)

        reader = threading.Thread(
            target=read_output,
            name=f"oneflow-output-{case.case_id}",
            daemon=True,
        )
        reader.start()
        reader_finished = False

        while process.poll() is None or not reader_finished:
            if cancel_event.is_set() and process.poll() is None:
                log(f"[{case.case_id}] Stop requested; terminating process.")
                process.terminate()
                try:
                    process.wait(timeout=3.0)
                except subprocess.TimeoutExpired:
                    process.kill()
                    process.wait()

            try:
                line = output_queue.get(timeout=0.1)
                if line is None:
                    reader_finished = True
                else:
                    log(f"[{case.case_id}] {line}")
            except queue.Empty:
                continue

        reader.join(timeout=1.0)
        return_code = process.returncode
        if cancel_event.is_set():
            state = RunState.CANCELLED
            message = "Run cancelled."
        elif return_code == 0:
            state = RunState.SUCCEEDED
            message = "Process exited successfully."
        else:
            state = RunState.FAILED
            message = f"OneFLOW exited with code {return_code}."

        return CaseRunResult(
            case.case_id, state, return_code, case_workspace, message
        )
    except (OSError, shutil.Error, ValueError) as error:
        return CaseRunResult(
            case.case_id, RunState.FAILED, None, workspace_root, str(error)
        )
