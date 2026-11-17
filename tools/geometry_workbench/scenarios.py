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



def save_scenario_contract(
    path: Path,
    scenario: Dict[str, Any],
    input_tables: Dict[str, Any],
    expected_tables: Dict[str, Any],
) -> None:
    """Persist edited contract tables while preserving scenario provenance.

    Execution results are deliberately not changed by editing or saving a
    contract. Only a real runner may record a result.
    """
    if not isinstance(input_tables, dict):
        raise ValueError("input tables must be a JSON object")
    if not isinstance(expected_tables, dict):
        raise ValueError("expected tables must be a JSON object")
    updated = dict(scenario)
    updated["input"] = input_tables
    updated["expected"] = expected_tables
    updated["result"] = scenario.get("result")
    path.write_text(
        json.dumps(updated, indent=2, ensure_ascii=True) + "\n",
        encoding="utf-8",
    )
