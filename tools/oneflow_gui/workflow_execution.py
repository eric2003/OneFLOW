"""Command construction for one-process, multi-directory OneFLOW workflows."""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path


def build_multi_directory_command(
    executable: Path,
    directories: Sequence[Path | str],
) -> list[str]:
    """Build a process argument list and preserve directory order.

    Paths are resolved so the command is independent of the process working
    directory. No shell interpolation is used.
    """
    if not directories:
        raise ValueError("A directory workflow must contain at least one directory.")

    resolved_directories: list[str] = []
    for index, directory in enumerate(directories):
        text = str(directory).strip()
        if not text:
            raise ValueError(f"Directory argument {index + 1} must not be empty.")
        resolved_directories.append(str(Path(text).expanduser().resolve()))

    return [str(Path(executable).expanduser().resolve()), "r", *resolved_directories]
