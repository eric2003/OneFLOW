"""Project-level operations for mesh layout templates."""

from __future__ import annotations

import shutil
from pathlib import Path


TEMPLATE_PROJECTS = {
    "rectangle": "minimal_rectangle",
    "flat_plate": "turbulent_flat_plate",
}


def ensure_template_project(
    repository_root: Path,
    projects_root: Path,
    kind: str,
) -> Path:
    """Copy a complete template project once, then preserve user-owned edits."""
    try:
        folder_name = TEMPLATE_PROJECTS[kind]
    except KeyError as error:
        raise ValueError(f"Unsupported mesh layout project type: {kind!r}") from error

    source = (repository_root / "grid" / folder_name).resolve()
    destination = (projects_root / folder_name).expanduser().resolve()
    if destination == source:
        raise ValueError("The editable project must not be the source template.")

    if destination.exists():
        if not destination.is_dir():
            raise FileExistsError(f"Project path exists and is not a directory: {destination}")
        expected_spec = "rectangle.json" if kind == "rectangle" else "flat_plate.json"
        if not (destination / expected_spec).is_file():
            raise FileExistsError(
                f"Project directory already exists but is not a complete template project: {destination}"
            )
        return destination

    if not source.is_dir():
        raise FileNotFoundError(f"Mesh layout template project not found: {source}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copytree(source, destination)
    return destination


__all__ = ["TEMPLATE_PROJECTS", "ensure_template_project"]
