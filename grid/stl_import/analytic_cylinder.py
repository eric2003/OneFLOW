"""Analytic lateral-cylinder construction between two user-defined axial planes."""
from __future__ import annotations
import argparse
import json
import math
import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from grid.stl_import.reader import Point, Triangle
from grid.stl_import.validation import analyze_surface


def build_cylinder_lateral_surface(radius: float, axis_start: float, axis_end: float,
                                   cut_start: float, cut_end: float,
                                   circumferential_divisions: int = 48,
                                   axial_divisions: int = 12) -> list[Triangle]:
    """Build only the analytic cylinder side between x=cut_start and x=cut_end."""
    values = (radius, axis_start, axis_end, cut_start, cut_end)
    if not all(math.isfinite(value) for value in values):
        raise ValueError("cylinder dimensions and cutting planes must be finite")
    if radius <= 0 or axis_end <= axis_start:
        raise ValueError("radius must be positive and axis_end must exceed axis_start")
    if not axis_start < cut_start < cut_end < axis_end:
        raise ValueError("cut planes must lie strictly inside the cylinder and be ordered")
    if circumferential_divisions < 3 or axial_divisions < 1:
        raise ValueError("use at least 3 circumferential and 1 axial division")
    rings: list[list[Point]] = []
    for i in range(axial_divisions + 1):
        x = cut_start + (cut_end - cut_start) * i / axial_divisions
        rings.append([
            (x, radius * math.cos(2.0 * math.pi * j / circumferential_divisions),
             radius * math.sin(2.0 * math.pi * j / circumferential_divisions))
            for j in range(circumferential_divisions)
        ])
    triangles: list[Triangle] = []
    for i in range(axial_divisions):
        for j in range(circumferential_divisions):
            k = (j + 1) % circumferential_divisions
            p00, p10 = rings[i][j], rings[i + 1][j]
            p01, p11 = rings[i][k], rings[i + 1][k]
            triangles.append((p00, p01, p10))
            triangles.append((p10, p01, p11))
    return triangles


def write_ascii_stl(path: Path, triangles: list[Triangle]) -> None:
    """Write a deterministic ASCII STL."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["solid analytic_cylinder_cut"]
    for a, b, c in triangles:
        ab = tuple(b[i] - a[i] for i in range(3))
        ac = tuple(c[i] - a[i] for i in range(3))
        normal = (ab[1]*ac[2]-ab[2]*ac[1], ab[2]*ac[0]-ab[0]*ac[2],
                  ab[0]*ac[1]-ab[1]*ac[0])
        magnitude = math.sqrt(sum(value * value for value in normal))
        if magnitude == 0:
            raise ValueError("cannot write a degenerate triangle")
        unit = tuple(value / magnitude for value in normal)
        lines.append("  facet normal " + " ".join(f"{v:.17g}" for v in unit))
        lines.append("    outer loop")
        for point in (a, b, c):
            lines.append("      vertex " + " ".join(f"{v:.17g}" for v in point))
        lines.extend(("    endloop", "  endfacet"))
    lines.append("endsolid analytic_cylinder_cut")
    path.write_text("\n".join(lines) + "\n", encoding="ascii", newline="\n")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-stl", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--radius", type=float, default=1.0)
    parser.add_argument("--axis-start", type=float, default=-2.0)
    parser.add_argument("--axis-end", type=float, default=2.0)
    parser.add_argument("--cut-start", type=float, required=True)
    parser.add_argument("--cut-end", type=float, required=True)
    parser.add_argument("--circumferential-divisions", type=int, default=48)
    parser.add_argument("--axial-divisions", type=int, default=12)
    args = parser.parse_args(argv)
    try:
        triangles = build_cylinder_lateral_surface(
            args.radius, args.axis_start, args.axis_end, args.cut_start, args.cut_end,
            args.circumferential_divisions, args.axial_divisions)
        write_ascii_stl(args.output_stl, triangles)
        quality = analyze_surface(triangles)
        report = {
            "operation": "analytic_cylinder_lateral_cut",
            "source_kind": "explicit_analytic_definition", "axis": "X",
            "center_yz": [0.0, 0.0], "radius": args.radius,
            "axis_extent": [args.axis_start, args.axis_end],
            "cut_planes": [{"equation": f"x = {args.cut_start:.17g}", "x": args.cut_start},
                           {"equation": f"x = {args.cut_end:.17g}", "x": args.cut_end}],
            "circumferential_divisions": args.circumferential_divisions,
            "axial_divisions": args.axial_divisions, "output_stl": str(args.output_stl),
            "surface_quality": quality, "end_caps_generated": False,
        }
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(json.dumps(report, indent=2) + "\n",
                               encoding="utf-8", newline="\n")
    except (OSError, ValueError) as exc:
        parser.exit(2, f"analytic_cylinder_cut: {exc}\n")
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
