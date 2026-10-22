"""Discover conventional OneFLOW executable locations without scanning the disk."""

from __future__ import annotations

import os
from pathlib import Path
from typing import Mapping


def executable_candidates(
    repository_root: Path,
    environ: Mapping[str, str] | None = None,
) -> list[Path]:
    """Return existing executable candidates in a predictable preference order.

    Repository build outputs are preferred for development. Conventional
    installation locations are fallbacks, and no recursive filesystem search
    is performed.
    """
    environment = os.environ if environ is None else environ
    root = repository_root.expanduser().resolve()

    candidates = [
        root / "build" / "bin" / "Debug" / "OneFLOW.exe",
        root / "build" / "bin" / "Release" / "OneFLOW.exe",
        Path("C:/dev/OneFLOW/bin/OneFLOW.exe"),
    ]

    program_files = environment.get("ProgramFiles")
    if program_files:
        candidates.append(Path(program_files) / "OneFLOW" / "bin" / "OneFLOW.exe")

    program_files_x86 = environment.get("ProgramFiles(x86)")
    if program_files_x86:
        candidates.append(
            Path(program_files_x86) / "OneFLOW" / "bin" / "OneFLOW.exe"
        )

    local_app_data = environment.get("LOCALAPPDATA")
    if local_app_data:
        candidates.append(
            Path(local_app_data) / "Programs" / "OneFLOW" / "bin" / "OneFLOW.exe"
        )

    # Always expose the repository's Debug and Release slots in the UI, even
    # before either configuration has been built. Installed fallbacks must exist.
    build_candidates = candidates[:2]
    fallback_candidates = candidates[2:]
    unique_candidates: list[Path] = []
    seen: set[str] = set()
    for candidate in build_candidates:
        normalized = str(candidate).casefold()
        if normalized not in seen:
            seen.add(normalized)
            unique_candidates.append(candidate.resolve())
    for candidate in fallback_candidates:
        normalized = str(candidate).casefold()
        if normalized in seen:
            continue
        seen.add(normalized)
        if candidate.is_file():
            unique_candidates.append(candidate.resolve())

    return unique_candidates
