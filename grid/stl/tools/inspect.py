"""Inspect an STL file before downstream grid conversion."""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from grid.stl.stl_to_gmsh import Stl2DError, extract_boundary, read_stl


def inspect(path: Path, tolerance: float = 1e-9) -> dict:
    triangles = read_stl(path)
    points = [point for triangle in triangles for point in triangle]
    bounds = {
        axis: [min(point[i] for point in points), max(point[i] for point in points)]
        for i, axis in enumerate(("x", "y", "z"))
    }
    z_values = bounds["z"]
    planar = z_values[1] - z_values[0] <= tolerance
    report = {
        "file": str(path),
        "format": "binary" if _is_binary_stl(path) else "ascii",
        "triangle_count": len(triangles),
        "coordinate_bounds": bounds,
        "planar_xy": planar,
        "constant_z": z_values[0] if planar else None,
    }
    try:
        loop = extract_boundary(triangles, tolerance)
    except Stl2DError as exc:
        report["topology_status"] = "unsupported_or_invalid"
        report["topology_diagnostic"] = str(exc)
    else:
        report["topology_status"] = "single_closed_planar_loop"
        report["boundary_vertex_count"] = len(loop)
        report["boundary_edge_count"] = len(loop)
    return report


def _is_binary_stl(path: Path) -> bool:
    data = path.read_bytes()
    if len(data) < 84:
        return False
    count = int.from_bytes(data[80:84], "little")
    return 84 + count * 50 == len(data)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    case_dir = Path(__file__).resolve().parents[1]
    parser.add_argument("--input", type=Path, default=case_dir / "input" / "body.stl")
    parser.add_argument("--output", type=Path, default=case_dir / "results" / "stl_report.json")
    parser.add_argument("--tolerance", type=float, default=1e-9)
    args = parser.parse_args(argv)
    try:
        report = inspect(args.input, args.tolerance)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8", newline="\n")
    except (OSError, Stl2DError, ValueError) as exc:
        parser.exit(2, f"stl_inspect: {exc}\n")
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
