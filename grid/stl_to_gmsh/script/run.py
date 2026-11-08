"""Run the standalone planar STL to Gmsh preprocessing case."""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

CASE_DIR = Path(__file__).resolve().parents[1]
REPOSITORY_ROOT = Path(__file__).resolve().parents[3]
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

from grid.stl.stl_to_gmsh import convert


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gmsh", help="optional Gmsh executable")
    args = parser.parse_args(argv)

    input_stl = CASE_DIR / "input" / "body.stl"
    config = CASE_DIR / "config" / "external_flow.json"
    results = CASE_DIR / "results"
    geo = results / "body.geo"
    msh = results / "body.msh"
    report_path = results / "mesh_report.json"

    results.mkdir(parents=True, exist_ok=True)
    report = convert(input_stl, config, geo)
    report["input"] = str(input_stl.relative_to(REPOSITORY_ROOT))
    report["physical_groups"] = [
        "Wall", "Inlet", "Outlet", "Symmetry", "Farfield", "Fluid"
    ]
    if args.gmsh:
        subprocess.run(
            [args.gmsh, str(geo), "-2", "-format", "msh2", "-o", str(msh)],
            check=True,
        )
        report["msh"] = str(msh.relative_to(REPOSITORY_ROOT))
    report_path.write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8", newline="\n"
    )
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
