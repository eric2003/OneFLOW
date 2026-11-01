"""Cross-platform process-group creation and process-tree termination."""

from __future__ import annotations

import os
import signal
import subprocess
from collections.abc import Callable
from typing import Any


CommandRunner = Callable[..., Any]


def process_group_options(*, platform_name: str | None = None) -> dict[str, object]:
    """Return Popen options that isolate a command and its descendants."""
    platform_value = (platform_name or os.name).lower()
    if platform_value in {"nt", "windows"}:
        return {"creationflags": getattr(subprocess, "CREATE_NEW_PROCESS_GROUP", 0x00000200)}
    return {"start_new_session": True}


def terminate_process_tree(
    process: subprocess.Popen[object],
    *,
    platform_name: str | None = None,
    command_runner: CommandRunner | None = None,
) -> bool:
    """Force-stop a process and its descendants; return whether it was running."""
    if process.poll() is not None:
        return False

    platform_value = (platform_name or os.name).lower()
    if platform_value in {"nt", "windows"}:
        runner = command_runner or subprocess.run
        try:
            result = runner(
                ["taskkill", "/PID", str(process.pid), "/T", "/F"],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
                check=False,
            )
        except OSError:
            process.terminate()
            return True
        if result.returncode != 0 and process.poll() is None:
            process.terminate()
        return True

    try:
        os.killpg(process.pid, signal.SIGKILL)
    except ProcessLookupError:
        pass
    except OSError:
        process.kill()
    return True
