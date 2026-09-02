#!/usr/bin/env python3
"""Validate the committed residual database against its source references."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


TEST_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(TEST_ROOT))

from residual_db import parse_residual_file, sha256_file, summarize  # noqa: E402


def verify(database_path: Path) -> list[str]:
    errors: list[str] = []
    database = json.loads(database_path.read_text(encoding="utf-8"))
    if database.get("schema") != "oneflow.residual-baseline":
        errors.append("unsupported or missing database schema")
    if database.get("schema_version") not in {1, 2, 3}:
        errors.append("unsupported or missing schema_version")

    suite_path = TEST_ROOT / database.get("suite", "")
    if not suite_path.is_file():
        errors.append(f"missing suite file: {suite_path}")
        return errors
    suite_cases = [
        line.strip()
        for line in suite_path.read_text(encoding="utf-8-sig").splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]
    database_cases = list(database.get("cases", {}))
    if database_cases != suite_cases:
        errors.append(
            f"database cases differ from suite: database={database_cases}, "
            f"suite={suite_cases}"
        )

    for case_name, case in database.get("cases", {}).items():
        if not case.get("files"):
            errors.append(f"{case_name}: no residual files")
        for relative_result, entry in case.get("files", {}).items():
            source = TEST_ROOT / entry.get("source_reference", "")
            label = f"{case_name}:{relative_result}"
            if not source.is_file():
                errors.append(f"{label}: missing source reference {source}")
                continue
            if sha256_file(source) != entry.get("sha256"):
                errors.append(f"{label}: source SHA256 differs")
            try:
                variables, rows = parse_residual_file(source)
            except (OSError, ValueError) as error:
                errors.append(f"{label}: {error}")
                continue
            if variables != entry.get("variables"):
                errors.append(f"{label}: variables differ")
            if rows != entry.get("rows"):
                errors.append(f"{label}: stored residual curve differs")
            if summarize(variables, rows) != entry.get("summary"):
                errors.append(f"{label}: summary differs")
    return errors


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "database",
        nargs="?",
        type=Path,
        default=Path(__file__).resolve().parent / "residual-baseline.json",
    )
    args = parser.parse_args()
    errors = verify(args.database.resolve())
    if errors:
        for error in errors:
            print(f"ERROR: {error}", file=sys.stderr)
        return 1
    database = json.loads(args.database.read_text(encoding="utf-8"))
    file_count = sum(
        len(case["files"]) for case in database.get("cases", {}).values()
    )
    print(
        f"residual database OK: {len(database.get('cases', {}))} cases, "
        f"{file_count} files, baseline={database.get('baseline_id')}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
