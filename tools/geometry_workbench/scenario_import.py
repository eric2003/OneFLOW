"""Bootstrap editable geometry scenarios from discovered GoogleTest records.

This creates traceable scenario skeletons; it deliberately does not pretend to
infer arbitrary C++ setup or assertions as executable JSON.
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any


def scenario_id_for_test(record: dict[str, Any]) -> str:
    if record.get("framework") != "GoogleTest":
        raise ValueError("only GoogleTest records can seed this scenario format")
    suite = str(record.get("suite") or "")
    name = str(record.get("name") or "")
    if not suite or not name:
        raise ValueError("GoogleTest record requires suite and name")
    slug = re.sub(r"[^a-z0-9]+", "_", f"{suite}_{name}".lower()).strip("_")
    return slug


def scenario_skeleton(record: dict[str, Any]) -> dict[str, Any]:
    """Build a provenance-preserving starter; never fabricate input or output."""
    scenario_id = scenario_id_for_test(record)
    return {
        "id": scenario_id,
        "title": f"{record['suite']}: {record['name']}",
        "description": "Generated from the repository GoogleTest inventory; complete the data contract from the test setup and assertions.",
        "test": {
            "framework": "GoogleTest",
            "suite": record["suite"],
            "name": record["name"],
            "source": record["source"],
        },
        "input": {},
        "expected": {},
        "result": None,
        "notes": [
            "This scenario skeleton was generated from GoogleTest source discovery.",
            "Input and expected values are intentionally empty: source discovery does not execute C++ or infer arbitrary assertions.",
            "Complete the JSON contract, then connect a real runner before recording result."
        ],
        "provenance": {
            "kind": "gtest-inventory-skeleton",
            "sourceLine": record.get("line"),
            "discovery": record.get("discovery", "source"),
        },
    }


def write_scenario_skeleton(
    record: dict[str, Any], scenario_dir: Path, *, overwrite: bool = False
) -> Path:
    scenario_dir.mkdir(parents=True, exist_ok=True)
    path = scenario_dir / f"{scenario_id_for_test(record)}.json"
    if path.exists() and not overwrite:
        raise FileExistsError(f"scenario already exists: {path}")
    path.write_text(
        json.dumps(scenario_skeleton(record), indent=2, ensure_ascii=True) + "\n",
        encoding="utf-8",
    )
    return path
