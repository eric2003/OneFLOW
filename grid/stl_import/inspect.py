"""Inspect STL geometry and topology without modifying the input."""
from __future__ import annotations
import argparse
import json
import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from grid.stl_import.reader import StlReadError, is_binary_stl, read_stl
from grid.stl_import.validation import analyze_surface


def inspect(path: Path, tolerance: float | None = None) -> dict:
    triangles = read_stl(path)
    surface = analyze_surface(triangles, tolerance)
    bounds = surface["coordinate_bounds"]
    z = bounds["z"]
    planar = z[1] - z[0] <= surface["weld_tolerance"]
    return {
        "file": str(path), "format": "binary" if is_binary_stl(path) else "ascii",
        "readable": True, "triangle_count": len(triangles),
        "coordinate_bounds": bounds, "planar_xy": planar,
        "constant_z": z[0] if planar else None, "surface_quality": surface,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--weld-tolerance", type=float, default=None,
                        help="absolute tolerance in model units; default is scale-aware")
    args = parser.parse_args(argv)
    try:
        report = inspect(args.input, args.weld_tolerance)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2) + "\n",
                               encoding="utf-8", newline="\n")
    except (OSError, StlReadError, ValueError) as exc:
        parser.exit(2, f"stl_inspect: {exc}\n")
    print(json.dumps(report, indent=2))
    return 0 if report["surface_quality"]["status"] != "invalid" else 1


if __name__ == "__main__":
    raise SystemExit(main())
