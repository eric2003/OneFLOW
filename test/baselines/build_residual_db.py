#!/usr/bin/env python3
"""Build the versioned residual baseline database from autotest outputs."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Dict


TEST_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(TEST_ROOT))

from residual_db import parse_residual_file, sha256_file, summarize  # noqa: E402


def build_database(suite_file: Path, validated_commit: str) -> Dict[str, Any]:
    cases: Dict[str, Any] = {}
    case_names = [
        line.strip()
        for line in suite_file.read_text(encoding="utf-8-sig").splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]
    for case_name in case_names:
        autotest_dir = TEST_ROOT / case_name / "autotest"
        files: Dict[str, Any] = {}
        for source in sorted(autotest_dir.glob("*.dat")):
            if source.name not in {"res.dat", "turbres.dat"}:
                continue
            variables, rows = parse_residual_file(source)
            files[f"results/{source.name}"] = {
                "source_reference": f"{case_name}/autotest/{source.name}",
                "sha256": sha256_file(source),
                "variables": variables,
                "rows": rows,
                "summary": summarize(variables, rows),
            }
        if not files:
            raise RuntimeError(f"{case_name}: no residual reference files found")
        cases[case_name] = {"files": files}

    return {
        "schema": "oneflow.residual-baseline",
        "schema_version": 1,
        "baseline_id": "cpu-serial-v1",
        "suite": str(suite_file.relative_to(TEST_ROOT)),
        "validated_commit": validated_commit,
        "platform": {
            "execution": "CPU serial",
            "launcher": "mpirun -np 1",
            "reference_outputs": "case/autotest/*.dat",
        },
        "comparison": {
            "index_columns": ["iter", "sub-iter"],
            "absolute_tolerance": 1.0e-8,
            "relative_tolerance": 0.0,
            "nan_policy": "reject",
        },
        "cases": cases,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--suite", type=Path, default=TEST_ROOT / "suites" / "cpu-serial.txt"
    )
    parser.add_argument("--validated-commit", default="working-tree")
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).resolve().parent / "residual-baseline.json",
    )
    args = parser.parse_args()
    database = build_database(args.suite.resolve(), args.validated_commit)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(database, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(f"wrote {args.output}")
    print(f"cases={len(database['cases'])}")
    print("residual_files=" + str(sum(
        len(case["files"]) for case in database["cases"].values()
    )))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
