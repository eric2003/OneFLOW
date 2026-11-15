"""Cross-platform filesystem path normalization and comparison helpers."""

from __future__ import annotations

import os
from pathlib import Path


def normalize_path(path: Path | str) -> Path:
    """Return an expanded, absolute path with existing symlinks resolved.

    Non-strict resolution also supports paths whose leaf or parent directories
    have not been created yet.
    """
    return Path(path).expanduser().resolve()


def paths_equal(left: Path | str, right: Path | str) -> bool:
    """Return whether two path spellings identify the same filesystem path.

    Existing paths are compared by filesystem identity first. This handles
    Windows 8.3 short names, case-insensitive aliases, and symlinks. The
    normalized, platform-aware lexical comparison is the fallback for paths
    that do not exist yet.
    """
    left_path = normalize_path(left)
    right_path = normalize_path(right)

    try:
        if os.path.samefile(left_path, right_path):
            return True
    except (OSError, ValueError, TypeError):
        pass

    left_key = os.path.normcase(os.path.normpath(os.fspath(left_path)))
    right_key = os.path.normcase(os.path.normpath(os.fspath(right_path)))
    return left_key == right_key


def is_path_within(path: Path | str, root: Path | str) -> bool:
    """Return whether path is equal to or below root.

    Both paths are resolved before the containment check so existing symlinks
    cannot be used to escape a protected workspace. Missing path components
    remain supported.
    """
    candidate = normalize_path(path)
    boundary = normalize_path(root)
    try:
        candidate.relative_to(boundary)
    except ValueError:
        return False
    return True


__all__ = ["is_path_within", "normalize_path", "paths_equal"]
