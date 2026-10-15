#!/usr/bin/env python3
"""Discover and validate registered OneFLOW grid cases."""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path
from typing import Any


INCLUDE_PATTERN = re.compile(r'"([^"]+)"\s*;')
CATALOG_SCHEMA = "oneflow.case-catalog"


def load_catalog(catalog_path: Path) -> dict[str, Any]:
    with catalog_path.open("r", encoding="utf-8-sig") as stream:
        catalog = json.load(stream)
    if catalog.get("schema") != CATALOG_SCHEMA:
        raise ValueError(f"Unsupported catalog schema: {catalog.get('schema')!r}")
    if catalog.get("version") != 1:
        raise ValueError(f"Unsupported catalog version: {catalog.get('version')!r}")
    cases = catalog.get("cases")
    if not isinstance(cases, list):
        raise ValueError("Catalog field 'cases' must be a list")

    seen_ids: set[str] = set()
    for entry in cases:
        if not isinstance(entry, dict):
            raise ValueError("Every catalog case must be an object")
        case_id = entry.get("id")
        directory = entry.get("directory")
        runtime_config = entry.get("runtimeConfig")
        if not isinstance(case_id, str) or not case_id:
            raise ValueError("Every catalog case must have a non-empty id")
        if case_id in seen_ids:
            raise ValueError(f"Duplicate case id: {case_id}")
        seen_ids.add(case_id)
        if not isinstance(directory, str) or not directory:
            raise ValueError(f"Case {case_id} must have a directory")
        if not isinstance(runtime_config, str) or not runtime_config:
            raise ValueError(f"Case {case_id} must have runtimeConfig")
    return catalog


def resolve_case_path(case_dir: Path, relative_path: str, label: str) -> Path:
    path = Path(relative_path)
    if path.is_absolute():
        raise ValueError(f"{label} must be relative to the case directory")
    resolved = (case_dir / path).resolve()
    try:
        resolved.relative_to(case_dir.resolve())
    except ValueError as error:
        raise ValueError(f"{label} escapes the case directory: {relative_path}") from error
    return resolved


def validate_case(repository_root: Path, entry: dict[str, Any]) -> list[str]:
    case_id = entry["id"]
    errors: list[str] = []
    case_dir = resolve_case_path(
        repository_root,
        entry["directory"],
        f"Case directory for {case_id}",
    )
    if not case_dir.is_dir():
        return [f"{case_id}: case directory does not exist: {case_dir}"]

    try:
        control_path = resolve_case_path(
            case_dir, entry["runtimeConfig"], f"runtimeConfig for {case_id}"
        )
        if not control_path.is_file():
            errors.append(f"{case_id}: runtime config not found: {control_path}")
        else:
            control_text = control_path.read_text(encoding="utf-8-sig")
            includes = INCLUDE_PATTERN.findall(control_text)
            if not includes:
                errors.append(f"{case_id}: runtime config contains no parameter includes")
            for include in includes:
                include_path = resolve_case_path(
                    control_path.parent, include, f"parameter include for {case_id}"
                )
                if not include_path.is_file():
                    errors.append(f"{case_id}: parameter file not found: {include_path}")

        layout = entry.get("layout")
        if layout:
            layout_path = resolve_case_path(case_dir, layout, f"layout for {case_id}")
            if not layout_path.is_file():
                errors.append(f"{case_id}: layout file not found: {layout_path}")

        profile = entry.get("profile")
        if profile:
            profile_path = resolve_case_path(case_dir, profile, f"profile for {case_id}")
            if not profile_path.is_file():
                errors.append(f"{case_id}: profile file not found: {profile_path}")
            else:
                with profile_path.open("r", encoding="utf-8-sig") as stream:
                    profile_data = json.load(stream)
                if profile_data.get("schema") != "oneflow.case":
                    errors.append(f"{case_id}: unsupported case profile schema")
                if profile_data.get("version") != 1:
                    errors.append(f"{case_id}: unsupported case profile version")
    except (OSError, ValueError, json.JSONDecodeError) as error:
        errors.append(f"{case_id}: {error}")
    return errors


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--root",
        type=Path,
        default=Path(__file__).resolve().parents[1],
        help="OneFLOW repository root",
    )
    parser.add_argument(
        "--catalog",
        type=Path,
        default=Path("grid/case_catalog.json"),
        help="catalog path relative to the repository root",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)
    subparsers.add_parser("list", help="list registered cases")
    show_parser = subparsers.add_parser("show", help="show one registered case")
    show_parser.add_argument("case_id")
    validate_parser = subparsers.add_parser("validate", help="validate catalog entries")
    validate_parser.add_argument("case_id", nargs="?", help="validate one case only")
    args = parser.parse_args()

    root = args.root.expanduser().resolve()
    catalog_path = args.catalog
    if not catalog_path.is_absolute():
        catalog_path = root / catalog_path
    try:
        catalog = load_catalog(catalog_path.resolve())
        entries = catalog["cases"]
        if args.command == "list":
            for entry in entries:
                print(f"{entry['id']}\t{entry.get('title', entry['id'])}")
            return 0

        if args.command == "show":
            matches = [entry for entry in entries if entry["id"] == args.case_id]
            if not matches:
                raise ValueError(f"Unknown case id: {args.case_id}")
            print(json.dumps(matches[0], indent=2, ensure_ascii=True))
            return 0

        selected = entries
        if args.case_id:
            selected = [entry for entry in entries if entry["id"] == args.case_id]
            if not selected:
                raise ValueError(f"Unknown case id: {args.case_id}")

        errors = [
            error
            for entry in selected
            for error in validate_case(root, entry)
        ]
        if errors:
            print("Case catalog validation failed:")
            for error in errors:
                print(f"  - {error}")
            return 1
        print(f"Case catalog validation passed: {len(selected)} case(s).")
        return 0
    except (OSError, ValueError, json.JSONDecodeError) as error:
        print(f"Mesh case command failed: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
