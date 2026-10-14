#!/usr/bin/env python3
"""Validate a OneFLOW case's control-file includes and grid-conversion paths."""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path


INCLUDE_PATTERN = re.compile(r'"([^"]+)"\s*;')
STRING_PATTERN = re.compile(
    r'^\s*string\s+([A-Za-z_]\w*)\s*=\s*"([^"]*)"\s*;?'
)
INTEGER_PATTERN = re.compile(r'^\s*int\s+([A-Za-z_]\w*)\s*=\s*(-?\d+)\b')


def strip_line_comment(line: str) -> str:
    """Remove C++-style comments without treating // inside strings as comments."""
    in_string = False
    escaped = False
    index = 0
    while index < len(line):
        char = line[index]
        if escaped:
            escaped = False
        elif char == "\\" and in_string:
            escaped = True
        elif char == '"':
            in_string = not in_string
        elif not in_string and line[index:index + 2] == "//":
            return line[:index]
        index += 1
    return line


def read_control(case_dir: Path) -> tuple[list[Path], list[str]]:
    control_path = case_dir / "script" / "control.txt"
    errors: list[str] = []
    if not control_path.is_file():
        return [], [f"Missing control file: {control_path}"]

    included: list[Path] = []
    text = control_path.read_text(encoding="utf-8-sig")
    for entry in INCLUDE_PATTERN.findall(text):
        path = case_dir / "script" / entry
        if not path.is_file():
            errors.append(f"Control include does not exist: script/{entry}")
        else:
            included.append(path)
    if not included and not errors:
        errors.append("script/control.txt contains no parameter-file includes")
    return included, errors


def read_parameters(files: list[Path]) -> tuple[dict[str, str], dict[str, int], list[str]]:
    strings: dict[str, str] = {}
    integers: dict[str, int] = {}
    errors: list[str] = []
    for path in files:
        for line_number, original in enumerate(
            path.read_text(encoding="utf-8-sig").splitlines(), start=1
        ):
            line = strip_line_comment(original).strip()
            match = STRING_PATTERN.match(line)
            if match:
                name, value = match.groups()
                if name in strings:
                    errors.append(
                        f"Duplicate string parameter '{name}' in {path.name}:{line_number}"
                    )
                strings[name] = value
                continue
            match = INTEGER_PATTERN.match(line)
            if match:
                name, value = match.groups()
                if name in integers:
                    errors.append(
                        f"Duplicate integer parameter '{name}' in {path.name}:{line_number}"
                    )
                integers[name] = int(value)
    return strings, integers, errors


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Check a OneFLOW case before running the executable."
    )
    parser.add_argument(
        "case_dir",
        type=Path,
        help="Case directory containing script/control.txt (repository template or workspace copy)",
    )
    parser.add_argument(
        "--create-output-dirs",
        action="store_true",
        help="Create missing parent directories for configured output files",
    )
    args = parser.parse_args()

    case_dir = args.case_dir.expanduser().resolve()
    if not case_dir.is_dir():
        print(f"ERROR: Case directory does not exist: {case_dir}", file=sys.stderr)
        return 2

    included, errors = read_control(case_dir)
    strings, integers, parameter_errors = read_parameters(included)
    errors.extend(parameter_errors)

    for name in ("simutask",):
        if name not in strings:
            errors.append(f"Missing required string parameter: {name}")
    if strings.get("simutask", "").lower() != "grid":
        errors.append("This validator currently supports simutask = Grid only")

    if "dimension" not in integers:
        errors.append("Missing required integer parameter: dimension")
    elif integers["dimension"] not in (1, 2, 3):
        errors.append("dimension must be 1, 2, or 3")

    if "gridObj" not in integers:
        errors.append("Missing required integer parameter: gridObj")
    elif integers["gridObj"] != 1:
        errors.append("This case validator currently supports gridObj = 1 (conversion) only")

    required_grid_strings = (
        "sourceGridType",
        "targetGridType",
        "gridLayoutFileName",
        "sourceGridFileName",
        "sourceGridBcName",
        "targetGridFileName",
    )
    for name in required_grid_strings:
        if not strings.get(name):
            errors.append(f"Missing required grid parameter: {name}")

    # GridConfig::FromDataBase reads gridLayoutFileName even for conversion;
    # the selected conversion workflow does not necessarily open that layout file.
    input_names = ("sourceGridFileName", "sourceGridBcName")
    for name in input_names:
        value = strings.get(name)
        if value:
            path = Path(value)
            resolved = path if path.is_absolute() else case_dir / path
            if not resolved.is_file():
                errors.append(f"{name} points to a missing input file: {resolved}")

    output_value = strings.get("targetGridFileName")
    if output_value:
        output_path = Path(output_value)
        output_path = output_path if output_path.is_absolute() else case_dir / output_path
        if args.create_output_dirs:
            output_path.parent.mkdir(parents=True, exist_ok=True)
        if not output_path.parent.is_dir():
            errors.append(
                f"Output directory does not exist for targetGridFileName: {output_path.parent} "
                "(use --create-output-dirs to create it)"
            )

    if errors:
        print("OneFLOW case validation failed:")
        for error in errors:
            print(f"  - {error}")
        return 1

    print("OneFLOW case validation passed.")
    print(f"  Case: {case_dir}")
    print(f"  Parameter files: {len(included)}")
    print(f"  Task: {strings['simutask']}")
    print(f"  Dimension: {integers['dimension']}")
    print(f"  Grid objective: {integers['gridObj']} (convert)")
    print(f"  Source: {strings['sourceGridType']} / {strings['sourceGridFileName']}")
    print(f"  Target: {strings['targetGridType']} / {strings['targetGridFileName']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
