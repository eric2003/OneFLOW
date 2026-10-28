"""Load geometry contract scenarios from JSON."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, List


def scenarios_dir(repository_root: Path | None = None) -> Path:
    if repository_root is None:
        return Path(__file__).resolve().parent / "scenarios"
    return (repository_root / "tools" / "geometry_workbench" / "scenarios").resolve()


def list_scenario_paths(repository_root: Path | None = None) -> List[Path]:
    root = scenarios_dir(repository_root)
    if not root.is_dir():
        return []
    return sorted(root.glob("*.json"))


def load_scenario(path: Path) -> Dict[str, Any]:
    with path.open("r", encoding="utf-8") as stream:
        return json.load(stream)
