#!/usr/bin/env python3
"""Versioned residual-baseline database reader and comparator."""

from __future__ import annotations

import hashlib
import json
import math
import re
from pathlib import Path
from typing import Any, Dict, List, Tuple


DEFAULT_DATABASE = (
    Path(__file__).resolve().parent / "baselines" / "residual-baseline.json"
)


def _variable_names(lines: List[str]) -> List[str]:
    names: List[str] = []
    in_variables = False
    for line in lines:
        stripped = line.strip()
        if stripped.startswith("Variables="):
            in_variables = True
            continue
        if in_variables:
            match = re.fullmatch(r'"([^"]+)"', stripped)
            if match:
                names.append(match.group(1))
                continue
            if stripped:
                break
    return names


def parse_residual_file(path: Path) -> Tuple[List[str], List[List[float]]]:
    """Parse a OneFLOW residual file into column names and numeric rows."""

    lines = path.read_text(encoding="utf-8-sig").splitlines()
    names = _variable_names(lines)
    rows: List[List[float]] = []
    expected_width = len(names)
    for line_number, line in enumerate(lines, start=1):
        fields = line.split()
        if not fields:
            continue
        try:
            row = [float(field) for field in fields]
        except ValueError:
            continue
        if expected_width and len(row) != expected_width:
            raise ValueError(
                f"{path}: line {line_number} has {len(row)} values; "
                f"expected {expected_width}"
            )
        rows.append(row)
    if not names:
        raise ValueError(f"{path}: missing Variables header")
    if not rows:
        raise ValueError(f"{path}: no numeric residual rows")
    return names, rows


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def summarize(names: List[str], rows: List[List[float]]) -> Dict[str, Any]:
    return {
        "row_count": len(rows),
        "first": dict(zip(names, rows[0])),
        "last": dict(zip(names, rows[-1])),
        "max_abs": {
            name: max(abs(row[index]) for row in rows)
            for index, name in enumerate(names[2:], start=2)
        },
    }


def load_database(path: Path | str = DEFAULT_DATABASE) -> Dict[str, Any]:
    with Path(path).open("r", encoding="utf-8") as stream:
        return json.load(stream)


def compare_case(
    case_directory: Path | str,
    database: Dict[str, Any] | Path | str = DEFAULT_DATABASE,
    absolute_tolerance: float | None = None,
    relative_tolerance: float | None = None,
) -> Dict[str, Any]:
    """Compare generated residual files with the case's JSON baseline."""

    case_path = Path(case_directory)
    db = load_database(database) if not isinstance(database, dict) else database
    case_name = case_path.name
    case_baseline = db.get("cases", {}).get(case_name)
    if case_baseline is None:
        return {"ok": False, "case": case_name, "errors": [
            f"no residual baseline entry for case '{case_name}'"
        ]}

    comparison = db.get("comparison", {})
    abs_tol = float(
        absolute_tolerance
        if absolute_tolerance is not None
        else comparison.get("absolute_tolerance", 1.0e-8)
    )
    rel_tol = float(
        relative_tolerance
        if relative_tolerance is not None
        else comparison.get("relative_tolerance", 0.0)
    )
    errors: List[str] = []
    files_checked = 0
    max_abs_diff = 0.0
    max_rel_diff = 0.0

    for relative_name, expected in case_baseline.get("files", {}).items():
        actual_path = case_path / relative_name
        if not actual_path.is_file():
            errors.append(f"missing residual file: {actual_path}")
            continue
        try:
            names, actual_rows = parse_residual_file(actual_path)
        except (OSError, ValueError) as error:
            errors.append(str(error))
            continue
        expected_names = expected["variables"]
        expected_rows = expected["rows"]
        if names != expected_names:
            errors.append(
                f"{actual_path}: variables differ: actual={names}, "
                f"expected={expected_names}"
            )
            continue
        if len(actual_rows) != len(expected_rows):
            errors.append(
                f"{actual_path}: row count differs: actual={len(actual_rows)}, "
                f"expected={len(expected_rows)}"
            )
            continue
        files_checked += 1
        for row_index, (actual_row, expected_row) in enumerate(
            zip(actual_rows, expected_rows), start=1
        ):
            for column_index, (actual, reference) in enumerate(
                zip(actual_row, expected_row)
            ):
                if column_index < 2:
                    if actual != reference:
                        errors.append(
                            f"{actual_path}: row {row_index}, index column "
                            f"{column_index} differs"
                        )
                    continue
                difference = abs(actual - reference)
                scale = max(abs(actual), abs(reference))
                tolerance = max(abs_tol, rel_tol * scale)
                max_abs_diff = max(max_abs_diff, difference)
                if scale:
                    max_rel_diff = max(max_rel_diff, difference / scale)
                if (
                    not math.isfinite(actual)
                    or not math.isfinite(reference)
                    or difference > tolerance
                ):
                    errors.append(
                        f"{actual_path}: row {row_index}, "
                        f"{names[column_index]} differs: actual={actual}, "
                        f"expected={reference}, tolerance={tolerance}"
                    )

    return {
        "ok": not errors,
        "case": case_name,
        "files_checked": files_checked,
        "max_absolute_difference": max_abs_diff,
        "max_relative_difference": max_rel_diff,
        "errors": errors,
    }
