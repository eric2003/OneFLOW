"""Qt-independent subprocess execution for workflow processes."""

from __future__ import annotations

import queue
import subprocess
import threading
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True, slots=True)
class ProcessResult:
    """Outcome of one external process."""

    status: str
    return_code: int
    message: str


def run_process(
    command: Sequence[str],
    *,
    cwd: Path,
    on_output: Callable[[str], None],
    is_cancelled: Callable[[], bool],
) -> ProcessResult:
    """Run a command, stream merged output, and support cooperative cancellation."""
    if not command or any(not isinstance(argument, str) or not argument for argument in command):
        raise ValueError("Command must contain non-empty string arguments.")

    process = subprocess.Popen(
        list(command),
        cwd=str(cwd),
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        errors="replace",
        bufsize=1,
    )
    output_queue: queue.Queue[str | None] = queue.Queue()

    def read_output() -> None:
        assert process.stdout is not None
        try:
            for line in process.stdout:
                output_queue.put(line.rstrip("\r\n"))
        finally:
            output_queue.put(None)

    reader = threading.Thread(target=read_output, daemon=True)
    reader.start()
    reader_finished = False
    cancellation_requested = False
    while process.poll() is None or not reader_finished:
        if is_cancelled() and process.poll() is None:
            cancellation_requested = True
            process.terminate()
            try:
                process.wait(timeout=3)
            except subprocess.TimeoutExpired:
                process.kill()
                process.wait()
        try:
            line = output_queue.get(timeout=0.1)
            if line is None:
                reader_finished = True
            else:
                on_output(line)
        except queue.Empty:
            continue
    reader.join(timeout=1)
    return_code = process.returncode
    if cancellation_requested:
        return ProcessResult("cancelled", return_code, "Workflow run cancelled.")
    if return_code == 0:
        return ProcessResult("succeeded", return_code, "OneFLOW exited successfully.")
    return ProcessResult("failed", return_code, f"OneFLOW exited with code {return_code}.")
