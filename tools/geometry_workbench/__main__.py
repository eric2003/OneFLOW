#!/usr/bin/env python3
"""CLI entry: python -m geometry_workbench from tools/ or via path."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from pipeline_stub import diff_tables, run_pipeline
from scenarios import load_scenario
from visualize import render_ascii


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Geometry contract workbench (table scenarios)")
    parser.add_argument(
        "--scenario",
        type=Path,
        default=Path(__file__).resolve().parent / "scenarios" / "two_zone_parent_split.json",
    )
    parser.add_argument("--backend", choices=("stub", "native"), default="stub")
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
        print("RESULT: PASS")
        code = 0
    orphans = expected.get("orphanParentGlobalIds") or []
    print()
    print(render_ascii(actual, orphan_ids=orphans))
    return code


if __name__ == "__main__":
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    raise SystemExit(main())
