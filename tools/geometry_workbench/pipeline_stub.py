"""Table-level interface-split pipeline: stub now, C++ binding later.

The stub hard-codes the golden result for known scenario ids so the
workbench UI can be developed before every production edge is exposed.
When RunSymmetricInterfaceSplitMappingPipeline is bound, replace
run_pipeline() body only — keep the same input/output dict shape.
"""

from __future__ import annotations

from copy import deepcopy
from typing import Any, Dict, List, Optional


def run_pipeline(scenario: Dict[str, Any], backend: str = "stub") -> Dict[str, Any]:
    """Run split mapping pipeline.

    backend:
      - "stub": return scenario["expected"] tables (development / UI track)
      - "native": reserved for future C++ binding; raises for now
    """
    if backend == "native":
        raise NotImplementedError(
            "native backend not wired yet; use backend=stub or bind "
            "RunSymmetricInterfaceSplitMappingPipeline"
        )
    if backend != "stub":
        raise ValueError(f"unknown backend: {backend}")

    expected = scenario.get("expected")
    if not expected:
        raise ValueError("scenario has no expected tables for stub backend")
    return deepcopy(expected)


def diff_tables(
    actual: Dict[str, Any],
    expected: Dict[str, Any],
    keys: Optional[List[str]] = None,
) -> List[str]:
    """Return human-readable mismatch lines; empty list means full match."""
    if keys is None:
        keys = sorted(set(expected) | set(actual))
    lines: List[str] = []
    for key in keys:
        if key not in expected:
            lines.append(f"{key}: unexpected key in actual")
            continue
        if key not in actual:
            lines.append(f"{key}: missing in actual")
            continue
        if actual[key] != expected[key]:
            lines.append(f"{key}: actual={actual[key]!r} expected={expected[key]!r}")
    return lines
