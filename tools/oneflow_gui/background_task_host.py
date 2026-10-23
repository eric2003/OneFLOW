"""Detached process host for a single OneFLOW command."""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from dataclasses import replace
from pathlib import Path
from typing import Callable

from .background_task_registry import (
    BackgroundTaskRegistry,
    utc_now,
)


def _write_log(log_file: Path, message: str) -> None:
    log_file.parent.mkdir(parents=True, exist_ok=True)
    with log_file.open("a", encoding="utf-8", newline="\n") as stream:
        stream.write(message.rstrip("\n") + "\n")
        stream.flush()


def run_host(
    *,
    root: Path,
    task_id: str,
    command: list[str],
    cwd: Path,
    stop_file: Path | None = None,
    poll_interval: float = 0.2,
    popen: Callable = subprocess.Popen,
) -> int:
    """Run a command independently of the GUI and persist its final state."""
    registry = BackgroundTaskRegistry(root)
    task = registry.load(task_id)
    log_file = Path(task.log_file)
    stop_path = stop_file or (registry.records_directory / f"{task_id}.stop")
    task = replace(
        task,
        host_pid=os.getpid(),
        process_pid=None,
        status="running",
        updated_at=utc_now(),
        message="Background host started.",
    )
    registry.save(task)
    _write_log(log_file, f"Background host PID: {os.getpid()}")
    process = None
    final_status = "failed"
    final_message = "Background host failed."
    return_code = 1
    try:
        with log_file.open("a", encoding="utf-8", newline="\n") as stream:
            stream.write("Command: " + json.dumps(command, ensure_ascii=True) + "\n")
            stream.flush()
            process = popen(
                command,
                cwd=str(cwd),
                stdin=subprocess.DEVNULL,
                stdout=stream,
                stderr=subprocess.STDOUT,
                shell=False,
                **({"creationflags": subprocess.CREATE_NEW_PROCESS_GROUP}
                   if os.name == "nt" else {"start_new_session": True}),
            )
            task = replace(
                task,
                process_pid=process.pid,
                updated_at=utc_now(),
                message="OneFLOW process started.",
            )
            registry.save(task)
            while process.poll() is None:
                if stop_path.exists():
                    task = replace(
                        task,
                        status="stopping",
                        updated_at=utc_now(),
                        message="Stop requested.",
                    )
                    registry.save(task)
                    try:
                        process.terminate()
                        process.wait(timeout=5)
                    except subprocess.TimeoutExpired:
                        process.kill()
                        process.wait()
                    break
                time.sleep(poll_interval)
            return_code = process.wait()
        if stop_path.exists():
            stop_path.unlink()
        if return_code == 0:
            final_status = "succeeded"
            final_message = "Process completed successfully."
        elif task.status == "stopping":
            final_status = "stopped"
            final_message = "Process stopped by request."
        else:
            final_status = "failed"
            final_message = f"Process exited with code {return_code}."
    except Exception as error:
        final_message = f"Background host error: {error}"
        _write_log(log_file, final_message)
        if process is not None and process.poll() is None:
            process.kill()
            process.wait()
    finally:
        current = registry.load(task_id)
        registry.save(
            replace(
                current,
                status=final_status,
                updated_at=utc_now(),
                message=final_message,
            )
        )
        _write_log(log_file, f"Final status: {final_status} ({final_message})")
    return return_code


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--task-id", required=True)
    parser.add_argument("--command-file", required=True, type=Path)
    parser.add_argument("--cwd", required=True, type=Path)
    args = parser.parse_args(argv)
    payload = json.loads(args.command_file.read_text(encoding="utf-8"))
    command = payload.get("command")
    if not isinstance(command, list) or not command or not all(
        isinstance(item, str) for item in command
    ):
        parser.error("command file must contain a non-empty string array named 'command'")
    return run_host(
        root=args.root,
        task_id=args.task_id,
        command=command,
        cwd=args.cwd,
    )


if __name__ == "__main__":
    raise SystemExit(main())
