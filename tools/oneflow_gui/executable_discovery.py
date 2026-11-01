"""Discover conventional OneFLOW executable locations without scanning the disk."""

from __future__ import annotations

import os
import shutil
from pathlib import Path
from typing import Mapping


def executable_candidates(
    repository_root: Path,
    environ: Mapping[str, str] | None = None,
    *,
    platform_name: str | None = None,
) -> list[Path]:
    """Return existing executable candidates in a predictable preference order.

    Repository build outputs are preferred for development. Conventional
    installation locations and PATH are fallbacks; no recursive search occurs.
    """
    environment = os.environ if environ is None else environ
    root = repository_root.expanduser().resolve()
    platform_value = (platform_name or os.name).lower()
    is_windows = platform_value in {"nt", "windows"}

    if is_windows:
        build_candidates = [
            root / "build" / "bin" / "Debug" / "OneFLOW.exe",
            root / "build" / "bin" / "Release" / "OneFLOW.exe",
        ]
        fallback_candidates: list[Path] = []
        for variable in ("ProgramFiles", "ProgramFiles(x86)"):
            directory = environment.get(variable)
            if directory:
                fallback_candidates.append(
                    Path(directory) / "OneFLOW" / "bin" / "OneFLOW.exe"
                )
        local_app_data = environment.get("LOCALAPPDATA")
        if local_app_data:
            fallback_candidates.append(
                Path(local_app_data) / "Programs" / "OneFLOW" / "bin" / "OneFLOW.exe"
            )
    else:
        build_candidates = [
            root / "build" / "bin" / "OneFLOW",
            root / "build" / "bin" / "Debug" / "OneFLOW",
            root / "build" / "bin" / "Release" / "OneFLOW",
        ]
        fallback_candidates = [
            Path("/usr/local/bin/OneFLOW"),
            Path("/usr/bin/OneFLOW"),
        ]
        path_executable = shutil.which("OneFLOW", path=environment.get("PATH"))
        if path_executable:
            fallback_candidates.append(Path(path_executable))

    unique_candidates: list[Path] = []
    seen: set[str] = set()
    for candidate in build_candidates:
        normalized = os.path.normcase(str(candidate))
        if normalized not in seen:
            seen.add(normalized)
            unique_candidates.append(candidate.resolve())
    for candidate in fallback_candidates:
        normalized = os.path.normcase(str(candidate))
        if normalized in seen:
            continue
        seen.add(normalized)
        if candidate.is_file():
            unique_candidates.append(candidate.resolve())

    return unique_candidates
