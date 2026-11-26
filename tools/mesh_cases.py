#!/usr/bin/env python3
"""Discover and validate registered OneFLOW grid cases."""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path
from typing import Any

try:
    from .case_lifecycle import validate_lifecycle_contract
except ImportError:  # Direct script execution from the repository root.
    from case_lifecycle import validate_lifecycle_contract


INCLUDE_PATTERN = re.compile(r'"([^"]+)"\s*;')
CATALOG_SCHEMA = "oneflow.case-catalog"
CASE_SCHEMA = "oneflow.case"


def load_catalog(catalog_path: Path) -> dict[str, Any]:
    with catalog_path.open("r", encoding="utf-8-sig") as stream:
        catalog = json.load(stream)
    if not isinstance(catalog, dict):
        raise ValueError("Catalog root must be an object")
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
        execution = entry.get("execution")
        python_case = isinstance(execution, dict) and execution.get("kind") == "python"
        if not isinstance(case_id, str) or not case_id:
            raise ValueError("Every catalog case must have a non-empty id")
        if case_id in seen_ids:
            raise ValueError(f"Duplicate case id: {case_id}")
        seen_ids.add(case_id)
        if not isinstance(directory, str) or not directory:
            raise ValueError(f"Case {case_id} must have a directory")
        dimension = entry.get("dimension")
        if dimension is not None and (not isinstance(dimension, int) or isinstance(dimension, bool) or dimension not in (2, 3)):
            raise ValueError(f"Case {case_id} dimension must be 2 or 3 when provided")
        if python_case:
            script = execution.get("script")
            if not isinstance(script, str) or not script:
                raise ValueError(f"Python case {case_id} must have execution.script")
        elif execution is not None:
            raise ValueError(f"Case {case_id} has an unsupported execution kind")
        elif runtime_config is not None and (
            not isinstance(runtime_config, str) or not runtime_config
        ):
            raise ValueError(
                f"Case {case_id} runtimeConfig must be a non-empty path when provided"
            )
        lifecycle_path_value = entry.get("lifecycle_contract")
        if lifecycle_path_value is not None:
            if not isinstance(lifecycle_path_value, str) or not lifecycle_path_value:
                raise ValueError(f"Case {case_id} lifecycle_contract must be a non-empty relative path")
            lifecycle_path = Path(lifecycle_path_value)
            if lifecycle_path.is_absolute() or ".." in lifecycle_path.parts:
                raise ValueError(f"Case {case_id} lifecycle_contract must stay within the grid directory")
            lifecycle_file = (catalog_path.parent / lifecycle_path).resolve()
            try:
                lifecycle_file.relative_to(catalog_path.parent.resolve())
            except ValueError as error:
                raise ValueError(f"Case {case_id} lifecycle_contract escapes the grid directory") from error
            if not lifecycle_file.is_file():
                raise ValueError(f"Case {case_id} lifecycle contract not found: {lifecycle_path_value}")
            with lifecycle_file.open("r", encoding="utf-8-sig") as stream:
                lifecycle = json.load(stream)
            lifecycle_errors = validate_lifecycle_contract(lifecycle, expected_case_id=case_id)
            if lifecycle_errors:
                details = "; ".join(lifecycle_errors)
                raise ValueError(f"Case {case_id} lifecycle contract is invalid: {details}")
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


def validate_profile(
    case_dir: Path,
    entry: dict[str, Any],
    profile_data: Any,
) -> list[str]:
    case_id = entry["id"]
    errors: list[str] = []
    if not isinstance(profile_data, dict):
        return [f"{case_id}: case profile root must be an object"]
    if profile_data.get("schema") != CASE_SCHEMA:
        errors.append(f"{case_id}: unsupported case profile schema")
    if profile_data.get("version") != 1:
        errors.append(f"{case_id}: unsupported case profile version")

    parameters = profile_data.get("parameters")
    if not isinstance(parameters, dict):
        errors.append(f"{case_id}: profile field 'parameters' must be an object")
        return errors
    runtime = parameters.get("runtime")
    if not isinstance(runtime, dict):
        errors.append(f"{case_id}: profile field 'parameters.runtime' must be an object")
        return errors

    profile_layout = runtime.get("layout")
    catalog_layout = entry.get("layout")
    if profile_layout != catalog_layout:
        errors.append(
            f"{case_id}: profile layout {profile_layout!r} does not match "
            f"catalog layout {catalog_layout!r}"
        )
    parameter_files = runtime.get("parameterFiles")
    if not isinstance(parameter_files, list) or not parameter_files:
        errors.append(f"{case_id}: runtime.parameterFiles must be a non-empty list")
    else:
        for parameter_file in parameter_files:
            if not isinstance(parameter_file, str) or not parameter_file:
                errors.append(f"{case_id}: runtime.parameterFiles entries must be non-empty strings")
                continue
            try:
                path = resolve_case_path(
                    case_dir, parameter_file, f"profile parameter file for {case_id}"
                )
                if not path.is_file():
                    errors.append(f"{case_id}: profile parameter file not found: {path}")
            except ValueError as error:
                errors.append(f"{case_id}: {error}")

    dimension = runtime.get("dimension")
    if not isinstance(dimension, int) or isinstance(dimension, bool) or dimension not in (1, 2, 3):
        errors.append(f"{case_id}: runtime.dimension must be 1, 2, or 3")

    grid_options = runtime.get("gridOptions")
    if grid_options is not None:
        if not isinstance(grid_options, dict):
            errors.append(f"{case_id}: runtime.gridOptions must be an object")
        else:
            if grid_options.get("topology") not in ("structured", "unstructured"):
                errors.append(
                    f"{case_id}: gridOptions.topology must be structured or unstructured"
                )
            if not isinstance(grid_options.get("multiBlock"), bool):
                errors.append(f"{case_id}: gridOptions.multiBlock must be a boolean")
            file_format = grid_options.get("fileFormat")
            if not isinstance(file_format, int) or isinstance(file_format, bool):
                errors.append(f"{case_id}: gridOptions.fileFormat must be an integer")
            scale = grid_options.get("scale")
            if not isinstance(scale, (int, float)) or isinstance(scale, bool):
                errors.append(f"{case_id}: gridOptions.scale must be numeric")
            translation = grid_options.get("translation")
            if (
                not isinstance(translation, list)
                or len(translation) != 3
                or any(not isinstance(value, (int, float)) or isinstance(value, bool)
                       for value in translation)
            ):
                errors.append(
                    f"{case_id}: gridOptions.translation must contain three numbers"
                )
            io_roots = grid_options.get("ioRoots")
            root_names = {
                "project", "input", "output", "meshInput", "flowInput",
                "meshOutput", "flowOutput",
            }
            if not isinstance(io_roots, dict):
                errors.append(f"{case_id}: gridOptions.ioRoots must be an object")
            else:
                for root_name in root_names:
                    if not isinstance(io_roots.get(root_name), str) or not io_roots[root_name]:
                        errors.append(
                            f"{case_id}: gridOptions.ioRoots.{root_name} must be a non-empty string"
                        )
    return errors


def validate_case(repository_root: Path, entry: dict[str, Any]) -> list[str]:
    case_id = entry["id"]
    errors: list[str] = []
    try:
        case_dir = resolve_case_path(
            repository_root,
            entry["directory"],
            f"Case directory for {case_id}",
        )
        if not case_dir.is_dir():
            return [f"{case_id}: case directory does not exist: {case_dir}"]

        execution = entry.get("execution")
        if execution is not None:
            if not isinstance(execution, dict) or execution.get("kind") != "python":
                errors.append(f"{case_id}: unsupported execution definition")
            else:
                script = execution.get("script")
                if not isinstance(script, str) or not script:
                    errors.append(f"{case_id}: Python execution requires a script")
                else:
                    script_path = resolve_case_path(
                        repository_root, script, f"Python execution script for {case_id}"
                    )
                    if not script_path.is_file():
                        errors.append(f"{case_id}: Python execution script not found: {script_path}")
                variants = execution.get("implementations", [])
                if not isinstance(variants, list):
                    errors.append(f"{case_id}: execution implementations must be an array")
                    variants = []
                variant_ids = []
                for variant in variants:
                    if not isinstance(variant, dict):
                        errors.append(f"{case_id}: execution implementation must be an object")
                        continue
                    variant_id = variant.get("id")
                    variant_kind = variant.get("kind")
                    variant_script = variant.get("script")
                    if not isinstance(variant_id, str) or not variant_id:
                        errors.append(f"{case_id}: execution implementation requires an id")
                    elif variant_id in variant_ids:
                        errors.append(f"{case_id}: duplicate implementation id {variant_id}")
                    else:
                        variant_ids.append(variant_id)
                    if variant_kind not in {"python", "cpp"}:
                        errors.append(f"{case_id}: unsupported implementation kind {variant_kind!r}")
                    if not isinstance(variant_script, str) or not variant_script:
                        errors.append(f"{case_id}: implementation {variant_id!r} requires a runner script")
                    else:
                        variant_path = resolve_case_path(
                            repository_root, variant_script,
                            f"runner script for {case_id}/{variant_id}",
                        )
                        if not variant_path.is_file():
                            errors.append(f"{case_id}: implementation runner not found: {variant_path}")
                default_implementation = execution.get("default_implementation")
                defaults = [item.get("id") for item in variants if isinstance(item, dict) and item.get("default")]
                if default_implementation is not None and default_implementation not in variant_ids:
                    errors.append(f"{case_id}: default implementation is not registered")
                if len(defaults) > 1:
                    errors.append(f"{case_id}: only one implementation may be marked default")
                if defaults and default_implementation and defaults[0] != default_implementation:
                    errors.append(f"{case_id}: conflicting default implementation declarations")
                required_files = ["expected.json"]
                if not execution.get("generates_input", False):
                    required_files.append("input/model.stl")
                for relative in required_files:
                    path = resolve_case_path(case_dir, relative, f"{relative} for {case_id}")
                    if not path.is_file():
                        errors.append(f"{case_id}: {relative} not found")
        else:
            runtime_config = entry.get("runtimeConfig")
            if runtime_config:
                control_path = resolve_case_path(
                    case_dir, runtime_config, f"runtimeConfig for {case_id}"
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
            else:
                script_dir = resolve_case_path(
                    case_dir, "script", f"script directory for {case_id}"
                )
                parameter_files = sorted(script_dir.glob("*.txt")) if script_dir.is_dir() else []
                if not parameter_files:
                    errors.append(f"{case_id}: no script/*.txt parameter files found")
                for parameter_file in parameter_files:
                    try:
                        parameter_file.resolve().relative_to(case_dir.resolve())
                    except ValueError:
                        errors.append(
                            f"{case_id}: parameter file escapes the case directory: {parameter_file}"
                        )
                    if not parameter_file.is_file():
                        errors.append(f"{case_id}: parameter file not found: {parameter_file}")

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
                errors.extend(validate_profile(case_dir, entry, profile_data))
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
