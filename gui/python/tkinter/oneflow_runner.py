"""Process execution primitives for the OneFLOW case runner.

This module deliberately contains no Tkinter code so command construction and
process orchestration can be reused by other interfaces and tested separately.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from threading import Event, Lock, Thread
from typing import Callable
import json
import subprocess


@dataclass(frozen=True)
class RunPlan:
    """An ordered list of independent OneFLOW process invocations."""

    executable: Path
    mode: str
    case_directories: tuple[Path, ...]
    stop_on_error: bool = True

    def validate(self) -> None:
        """Reject invalid plans before starting any external process."""
        if not self.executable.is_file():
            raise ValueError(f"OneFLOW executable does not exist: {self.executable}")
        if not self.mode.strip():
            raise ValueError("Execution mode cannot be empty.")
        if not self.case_directories:
            raise ValueError("Add at least one case directory before running.")

        for case_directory in self.case_directories:
            if not case_directory.is_dir():
                raise ValueError(f"Case directory does not exist: {case_directory}")

    def command_for(self, case_directory: Path) -> tuple[str, ...]:
        """Build argv without a shell, preserving spaces in Windows paths."""
        return (
            str(self.executable.resolve()),
            self.mode,
            str(case_directory.resolve()),
        )


def save_run_plan(path: Path, plan: RunPlan) -> None:
    """Write a portable, versioned JSON description of the current run order."""
    payload = {
        "schema": "oneflow.run-plan",
        "version": 1,
        "executable": str(plan.executable),
        "mode": plan.mode,
        "case_directories": [str(item) for item in plan.case_directories],
        "stop_on_error": plan.stop_on_error,
    }
    path.write_text(
        json.dumps(payload, indent=2, ensure_ascii=True) + "\n",
        encoding="utf-8",
    )


def load_run_plan(path: Path) -> RunPlan:
    """Read a saved plan and resolve relative paths against the plan file."""
    payload = json.loads(path.read_text(encoding="utf-8"))
    if payload.get("schema") != "oneflow.run-plan" or payload.get("version") != 1:
        raise ValueError(f"Unsupported run-plan format: {path}")

    base_directory = path.resolve().parent

    def resolve_saved_path(value: str) -> Path:
        candidate = Path(value).expanduser()
        if not candidate.is_absolute():
            candidate = base_directory / candidate
        return candidate.resolve()

    return RunPlan(
        executable=resolve_saved_path(payload["executable"]),
        mode=str(payload.get("mode", "r")),
        case_directories=tuple(
            resolve_saved_path(item) for item in payload.get("case_directories", [])
        ),
        stop_on_error=bool(payload.get("stop_on_error", True)),
    )


EventSink = Callable[[dict], None]


class OneFlowRunner:
    """Run selected cases sequentially and stream process output as events.

    The runner never calls GUI methods. Its event sink should enqueue events
    for the UI thread instead of updating widgets from this worker thread.
    """

    def __init__(self, event_sink: EventSink) -> None:
        self._event_sink = event_sink
        self._stop_requested = Event()
        self._process_lock = Lock()
        self._process: subprocess.Popen[str] | None = None
        self._thread: Thread | None = None

    @property
    def is_running(self) -> bool:
        return self._thread is not None and self._thread.is_alive()

    def start(self, plan: RunPlan) -> None:
        """Validate and start a run plan on a dedicated worker thread."""
        plan.validate()
        if self.is_running:
            raise RuntimeError("A OneFLOW run is already active.")

        self._stop_requested.clear()
        self._thread = Thread(
            target=self._run_plan,
            args=(plan,),
            name="oneflow-case-runner",
            daemon=True,
        )
        self._thread.start()

    def stop(self) -> None:
        """Request cancellation and terminate the currently running process."""
        self._stop_requested.set()
        with self._process_lock:
            process = self._process
        if process is not None and process.poll() is None:
            process.terminate()

    def _emit(self, event_type: str, **payload: object) -> None:
        self._event_sink({"type": event_type, **payload})

    def _run_plan(self, plan: RunPlan) -> None:
        final_code = 0
        stopped = False
        total_cases = len(plan.case_directories)

        try:
            for index, case_directory in enumerate(plan.case_directories):
                if self._stop_requested.is_set():
                    stopped = True
                    final_code = 130
                    break

                command = plan.command_for(case_directory)
                self._emit(
                    "case_started",
                    index=index,
                    total=total_cases,
                    case=str(case_directory),
                    command=list(command),
                )

                try:
                    process = subprocess.Popen(
                        command,
                        cwd=str(case_directory),
                        stdout=subprocess.PIPE,
                        stderr=subprocess.STDOUT,
                        text=True,
                        encoding="utf-8",
                        errors="replace",
                        bufsize=1,
                    )
                except OSError as error:
                    final_code = 1
                    self._emit(
                        "error",
                        case=str(case_directory),
                        message=f"Could not start OneFLOW: {error}",
                    )
                    break

                with self._process_lock:
                    self._process = process
                    stop_after_start = self._stop_requested.is_set()

                # Close the small race where Stop is pressed just after Popen
                # returns but before the process is published to the UI thread.
                if stop_after_start:
                    process.terminate()

                assert process.stdout is not None
                for line in process.stdout:
                    self._emit("log", case=str(case_directory), line=line)

                return_code = process.wait()
                with self._process_lock:
                    self._process = None

                self._emit(
                    "case_finished",
                    index=index,
                    total=total_cases,
                    case=str(case_directory),
                    return_code=return_code,
                )

                if self._stop_requested.is_set():
                    stopped = True
                    final_code = 130
                    break

                if return_code != 0:
                    final_code = return_code
                    if plan.stop_on_error:
                        break
        except Exception as error:  # Report worker failures to the UI thread.
            final_code = 1
            self._emit("error", message=f"Runner failed: {error}")
        finally:
            with self._process_lock:
                self._process = None
            self._emit(
                "finished",
                return_code=final_code,
                stopped=stopped,
            )
