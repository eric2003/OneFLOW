#!/usr/bin/env python3
"""CLI entry for geometry contract scenarios.

Default: text only (ASCII + PASS/FAIL). Optional static PNG via --plot-file.
No interactive window is opened unless explicitly requested later.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from pipeline_stub import diff_tables, run_pipeline
from scenarios import load_scenario
from visualize import render_ascii, save_topology_png


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Geometry contract workbench (table scenarios)",
        allow_abbrev=False,
    )
    parser.add_argument(
        "--scenario",
        type=Path,
        default=Path(__file__).resolve().parent / "scenarios" / "two_zone_parent_split.json",
    )
    parser.add_argument("--backend", choices=("stub", "native"), default="stub")
    parser.add_argument(
        "--plot-file",
        type=Path,
        default=None,
        help="Optional static PNG path (matplotlib Agg; no interactive window)",
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
        print("RESULT: PASS")
        code = 0
    orphans = expected.get("orphanParentGlobalIds") or []
    print()
    print(render_ascii(actual, orphan_ids=orphans))
    if args.plot_file is not None:
        try:
            out = save_topology_png(actual, args.plot_file, orphan_ids=orphans)
            print(f"\nwrote static figure: {out}")
        except ImportError:
            print("\n--plot-file ignored: matplotlib not installed", file=sys.stderr)
            code = code or 0
    return code


if __name__ == "__main__":
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    raise SystemExit(main())
