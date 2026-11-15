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


def _existing_ancestor(path: Path) -> tuple[Path, tuple[str, ...]]:
    """Split a normalized path into its nearest existing ancestor and suffix."""
    current = path
    suffix: list[str] = []
    while not current.exists():
        parent = current.parent
        if parent == current:
            break
        suffix.append(current.name)
        current = parent
    return current, tuple(reversed(suffix))


def _same_file(left: Path, right: Path) -> bool:
    try:
        return os.path.samefile(left, right)
    except (OSError, ValueError, TypeError):
        return False


def paths_equal(left: Path | str, right: Path | str) -> bool:
    """Return whether two path spellings identify the same filesystem path.

    Existing paths are compared by filesystem identity first. For missing
    descendants, the nearest existing ancestors are compared by identity and
    the remaining path components are compared with platform-aware casing.
    This also handles Windows 8.3 aliases when the final target does not exist.
    """
    left_path = normalize_path(left)
    right_path = normalize_path(right)

    if _same_file(left_path, right_path):
        return True

    left_ancestor, left_suffix = _existing_ancestor(left_path)
    right_ancestor, right_suffix = _existing_ancestor(right_path)
    if _same_file(left_ancestor, right_ancestor):
        left_parts = tuple(os.path.normcase(part) for part in left_suffix)
        right_parts = tuple(os.path.normcase(part) for part in right_suffix)
        if left_parts == right_parts:
            return True

    left_key = os.path.normcase(os.path.normpath(os.fspath(left_path)))
    right_key = os.path.normcase(os.path.normpath(os.fspath(right_path)))
    return left_key == right_key


def is_path_within(path: Path | str, root: Path | str) -> bool:
    """Return whether path is equal to or below root.

    Both paths are resolved before the containment check so existing symlinks
    cannot be used to escape a protected workspace. Ancestor identity checks
    also support missing descendants reached through path aliases.
    """
    candidate = normalize_path(path)
    boundary = normalize_path(root)
    current = candidate
    while True:
        if paths_equal(current, boundary):
            return True
        parent = current.parent
        if parent == current:
            return False
        current = parent


__all__ = ["is_path_within", "normalize_path", "paths_equal"]
