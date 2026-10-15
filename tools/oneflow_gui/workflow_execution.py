"""Command construction for one-process, multi-directory OneFLOW workflows."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
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



def read_multi_directory_step(steps: Sequence[Mapping[str, object]]) -> list[str]:
    """Read and validate the supported multi-directory workflow step."""
    matching_steps = [
        step for step in steps if step.get("type") == "multi_directory"
    ]
    if len(matching_steps) != 1:
        raise ValueError(
            "A workflow must contain exactly one multi_directory step."
        )

    directories = matching_steps[0].get("directories")
    if not isinstance(directories, (list, tuple)) or not directories:
        raise ValueError(
            "A multi_directory step must contain at least one directory."
        )

    normalized: list[str] = []
    for index, directory in enumerate(directories, start=1):
        if not isinstance(directory, str) or not directory.strip():
            raise ValueError(
                f"Directory argument {index} must be a non-empty string."
            )
        normalized.append(directory.strip())
    return normalized
