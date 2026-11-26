"""Run the smallest safe Python test slice for the current repository changes.

The gate reuses the existing Python execution registry. Unknown source areas
fall back to the complete registered Python suite rather than being skipped.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path
from typing import Iterable

from tools.oneflow_gui.test_runner import _PYTHON_UNIT_TESTS
from .preflight import _changed_files, collect_preflight_issues, format_report


def _module(arguments: tuple[str, ...]) -> str | None:
    if len(arguments) >= 3 and arguments[:2] == ("-m", "unittest"):
        if arguments[2] == "discover":
            return None
        return arguments[2]
    return None


def select_test_commands(
    changed_files: Iterable[str],
    registry: tuple = _PYTHON_UNIT_TESTS,
) -> tuple[tuple[str, tuple[str, ...]], ...]:
    """Select registered test commands; unknown code paths trigger full suite."""
    paths = tuple(sorted({str(path).replace("\\", "/").lstrip("./") for path in changed_files}))
    if not paths:
        return ()

    scopes: set[str] = set()
    full_suite = False
    for path in paths:
        if path.startswith(("grid/stl_to_gmsh/", "grid/tests/python/stl_to_gmsh/")):
            scopes.update(("grid.tests.python.stl_to_gmsh.", "grid.tests.python.stl_import."))
            scopes.update((
                "tools.oneflow_gui.test_mesh_validation",
                "tools.oneflow_gui.test_geometry_asset_store",
                "tools.oneflow_gui.test_mesh_generation_dialog",
            ))
        elif path.startswith(("grid/stl_import/", "grid/tests/python/stl_import/")):
            scopes.update(("grid.tests.python.stl_import.", "grid.tests.python.stl_to_gmsh."))
            scopes.update(("tools.oneflow_gui.test_cylinder_cut_dialog",
                           "tools.oneflow_gui.test_oblique_cylinder_cut_dialog"))
        elif path.startswith("tools/verification/"):
            scopes.add("tools.verification.")
        elif path.startswith("tools/oneflow_gui/"):
            scopes.add("tools.oneflow_gui.")
        elif path.startswith(("tools/mesh_layout/", "grid/tests/python/mesh_layout/")):
            scopes.add("grid.tests.python.mesh_layout.")
            scopes.update(("tools.oneflow_gui.test_mesh_layout_projects",
                           "tools.oneflow_gui.test_mesh_generation_dialog"))
        elif path.startswith(("grid/", "tools/", "tests/", "codes/")):
            full_suite = True
        elif path.startswith((".github/", "docs/", "README", "AGENTS.md")):
            continue
        else:
            # Unknown non-documentation paths may affect Python services.
            full_suite = True

    selected = []
    for label, arguments in registry:
        module = _module(tuple(arguments))
        if full_suite or (module and any(
            module == scope or module.startswith(scope) for scope in scopes
        )):
            selected.append((label, tuple(arguments)))
    return tuple(selected)


def run_change_gate(
    root: Path,
    changed_files: Iterable[str],
    *,
    dry_run: bool = False,
    registry: tuple = _PYTHON_UNIT_TESTS,
) -> int:
    """Run repository preflight and affected registered Python test commands."""
    root = Path(root).resolve()
    changed = tuple(sorted(set(changed_files)))
    issues = collect_preflight_issues(root)
    print(format_report(issues, changed), flush=True)
    if issues:
        return 1

    selected = select_test_commands(changed, registry)
    if not selected:
        print("Change Gate: no affected Python tests selected (no code changes).")
        return 0

    print(f"Change Gate: {len(selected)} registered Python test command(s) selected.")
    failures = 0
    for index, (label, arguments) in enumerate(selected, 1):
        command = [sys.executable, *arguments]
        rendered = subprocess.list2cmdline(command) if sys.platform == "win32" else " ".join(command)
        print(f"[{index}/{len(selected)}] {label}\n$ {rendered}", flush=True)
        if dry_run:
            continue
        result = subprocess.run(command, cwd=root, check=False)
        if result.returncode:
            failures += 1
            print(f"FAIL ({result.returncode}): {label}", flush=True)
        else:
            print(f"PASS: {label}", flush=True)

    if dry_run:
        print("DRY RUN: commands were planned but not executed.")
        return 0
    print(f"Change Gate: {len(selected) - failures} passed, {failures} failed.")
    return 1 if failures else 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Run repository consistency checks and affected registered Python tests."
    )
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--changed", action="store_true",
                        help="detect staged, unstaged, and untracked Git paths")
    parser.add_argument("--all", action="store_true",
                        help="run the complete registered Python test suite")
    parser.add_argument("--dry-run", action="store_true",
                        help="print selected commands without executing them")
    args = parser.parse_args(argv)
    if args.changed == args.all:
        parser.error("choose exactly one of --changed or --all")
    root = args.root.resolve()
    try:
        changed = _changed_files(root) if args.changed else ("tools/change_gate_full_suite",)
        return run_change_gate(root, changed, dry_run=args.dry_run)
    except (OSError, RuntimeError) as exc:
        print(f"Change Gate: FAIL: {exc}")
        return 2


if __name__ == "__main__":
    sys.exit(main())
