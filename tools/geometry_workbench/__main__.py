#!/usr/bin/env python3
"""Geometry workbench CLI: edit-oriented table scenarios, not gtest launcher.

Examples:
  python -m geometry_workbench --scenario scenarios/two_zone_parent_split.json
  python -m geometry_workbench --scenario ... --plot plotly
  python -m geometry_workbench --scenario ... --plot none
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from pipeline_stub import diff_tables, run_pipeline
from visualize import render


def load_scenario(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as stream:
        return json.load(stream)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Table-level geometry workbench (stub pipeline)")
    parser.add_argument(
        "--scenario",
        type=Path,
        default=Path(__file__).resolve().parent / "scenarios" / "two_zone_parent_split.json",
        help="JSON scenario with input + expected tables",
    )
    parser.add_argument(
        "--backend",
        choices=("stub", "native"),
        default="stub",
        help="stub uses expected tables; native reserved for C++ binding",
    )
    parser.add_argument(
        "--plot",
        choices=("auto", "plotly", "matplotlib", "none"),
        default="auto",
        help="visualization backend (auto tries plotly then matplotlib)",
    )
    args = parser.parse_args(argv)

    scenario = load_scenario(args.scenario)
    actual = run_pipeline(scenario, backend=args.backend)
    expected = scenario.get("expected") or {}
    mismatches = diff_tables(actual, expected)

    print(f"scenario: {scenario.get('id', args.scenario.name)}")
    print(f"backend:  {args.backend}")
    if mismatches:
        print("RESULT: FAIL")
        for line in mismatches:
            print(f"  - {line}")
        code = 1
    else:
        print("RESULT: PASS (actual tables match expected)")
        code = 0

    orphans = expected.get("orphanParentGlobalIds") or []
    print()
    print(render(actual, orphan_ids=orphans, prefer=args.plot))
    return code


if __name__ == "__main__":
    # Allow `python tools/geometry_workbench/__main__.py` from repo root
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    raise SystemExit(main())
