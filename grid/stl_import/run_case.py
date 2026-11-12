"""Run one catalogued STL reader compatibility case."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from grid.stl_import.reader import StlReadError, is_binary_stl, read_stl
from grid.stl_import.surface import StlGeometryError, extract_boundary


def _format(path: Path) -> str:
    data = path.read_bytes()
    if is_binary_stl(path):
        return "binary"
    try:
        data.decode("ascii")
    except UnicodeDecodeError:
        return "unknown"
    return "ascii"


def run_case(case_directory: Path) -> dict:
    case_directory = case_directory.resolve()
    source = case_directory / "input" / "model.stl"
    expected = json.loads((case_directory / "expected.json").read_text(encoding="utf-8"))
    report = {"case": case_directory.name, "input": "input/model.stl", "format": _format(source),
              "readable": False, "triangle_count": None, "topology": "not_checked"}
    try:
        triangles = read_stl(source)
    except (OSError, StlReadError) as error:
        report["read_error"] = str(error)
    else:
        report["readable"] = True
        report["triangle_count"] = len(triangles)
        try:
            boundary = extract_boundary(triangles)
        except StlGeometryError as error:
            report["topology"] = "invalid"
            report["topology_error"] = str(error)
        else:
            report["topology"] = "valid"
            report["boundary_vertex_count"] = len(boundary)

    mismatches = []
    for key in ("readable", "format", "triangle_count", "topology"):
        if key in expected and report.get(key) != expected[key]:
            mismatches.append(f"{key}: expected {expected[key]!r}, got {report.get(key)!r}")
    if "read_error_contains" in expected and expected["read_error_contains"] not in report.get("read_error", ""):
        mismatches.append("read_error does not contain the expected diagnostic")
    if "topology_error_contains" in expected and expected["topology_error_contains"] not in report.get("topology_error", ""):
        mismatches.append("topology_error does not contain the expected diagnostic")
    report["passed"] = not mismatches
    report["expectation_mismatches"] = mismatches
    output = case_directory / "results" / "stl_report.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8", newline="\n")
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case-directory", type=Path, required=True)
    args = parser.parse_args()
    try:
        report = run_case(args.case_directory)
    except (OSError, ValueError, json.JSONDecodeError) as error:
        print(f"STL case runner failed: {error}", file=sys.stderr)
        return 2
    print(json.dumps(report, indent=2))
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
