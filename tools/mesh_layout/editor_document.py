"""Qt-independent document operations for the mesh layout editor."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any


def load_specification(path: Path) -> dict[str, Any]:
    """Load and validate a supported mesh layout specification."""
    spec = json.loads(path.read_text(encoding="utf-8"))
    compile_specification(spec)
    return spec


def compile_specification(spec: dict[str, Any]) -> tuple[str, str]:
    """Return the legacy layout text and SVG preview for a supported specification."""
    kind = spec.get("type")
    if kind == "rectangle":
        from tools.mesh_layout.compile_rectangle import (
            compile_rectangle,
            render_rectangle_svg,
        )

        return compile_rectangle(spec), render_rectangle_svg(spec)
    if kind == "flat_plate":
        from tools.mesh_layout.compile_flat_plate import (
            compile_flat_plate,
            render_flat_plate_svg,
        )

        return compile_flat_plate(spec), render_flat_plate_svg(spec)
    raise ValueError(f"Unsupported mesh layout type: {kind!r}")


def save_specification(path: Path, spec: dict[str, Any]) -> None:
    """Validate and save a specification as readable UTF-8 JSON."""
    compile_specification(spec)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(spec, indent=2, ensure_ascii=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )


__all__ = ["compile_specification", "load_specification", "save_specification"]
