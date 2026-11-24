"""Command-line adapter for blind-input generic STL surface cutting."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from grid.stl_import.reader import read_stl
from grid.stl_import.surface_cut import cut_closed_surface_between_planes


def _vector(values: list[float], label: str) -> tuple[float, float, float]:
    if len(values) != 3 or not all(math.isfinite(value) for value in values):
        raise ValueError(f"{label} requires three finite coordinates")
    return (values[0], values[1], values[2])


def _write_ascii_stl(path: Path, triangles, solid_name: str) -> None:
    with path.open("w", encoding="ascii", newline="\n") as stream:
        stream.write(f"solid {solid_name}\n")
        for a, b, c in triangles:
            ab = (b[0]-a[0], b[1]-a[1], b[2]-a[2])
            ac = (c[0]-a[0], c[1]-a[1], c[2]-a[2])
            normal = (
                ab[1]*ac[2]-ab[2]*ac[1],
                ab[2]*ac[0]-ab[0]*ac[2],
                ab[0]*ac[1]-ab[1]*ac[0],
            )
            length = math.sqrt(sum(value*value for value in normal))
            if length:
                normal = tuple(value/length for value in normal)
            stream.write("  facet normal {:.9g} {:.9g} {:.9g}\n".format(*normal))
            stream.write("    outer loop\n")
            for point in (a, b, c):
                stream.write("      vertex {:.9g} {:.9g} {:.9g}\n".format(*point))
            stream.write("    endloop\n  endfacet\n")
        stream.write(f"endsolid {solid_name}\n")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output-directory", required=True, type=Path)
    parser.add_argument("--inlet-point", nargs=3, required=True, type=float)
    parser.add_argument("--inlet-normal", nargs=3, required=True, type=float)
    parser.add_argument("--outlet-point", nargs=3, required=True, type=float)
    parser.add_argument("--outlet-normal", nargs=3, required=True, type=float)
    parser.add_argument("--inlet-keep", choices=("positive", "negative"), default="positive")
    parser.add_argument("--outlet-keep", choices=("positive", "negative"), default="negative")
    parser.add_argument("--tolerance", type=float, default=1.0e-8)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    source_path = args.input.expanduser().resolve()
    output = args.output_directory.expanduser().resolve()
    if not source_path.is_file():
        raise SystemExit(f"Input STL does not exist: {source_path}")
    output.mkdir(parents=True, exist_ok=True)
    source = read_stl(source_path)
    result = cut_closed_surface_between_planes(
        source,
        _vector(args.inlet_point, "inlet point"),
        _vector(args.inlet_normal, "inlet normal"),
        _vector(args.outlet_point, "outlet point"),
        _vector(args.outlet_normal, "outlet normal"),
        tolerance=args.tolerance,
        inlet_keep_positive=args.inlet_keep == "positive",
        outlet_keep_positive=args.outlet_keep == "positive",
    )
    artifacts = {
        "wall": result["side_triangles"],
        "inlet": result["inlet_cap_triangles"],
        "outlet": result["outlet_cap_triangles"],
        "closed_surface": result["closed_triangles"],
    }
    paths = {}
    for name, triangles in artifacts.items():
        path = output / f"{name}.stl"
        _write_ascii_stl(path, triangles, name)
        paths[name] = str(path)
    import hashlib
    region_manifest = {
        "schema": "oneflow.stl-surface-regions",
        "version": 1,
        "closed_surface_file": "closed_surface.stl",
        "closed_surface_sha256": hashlib.sha256((output / "closed_surface.stl").read_bytes()).hexdigest(),
        "regions": {
            name: {
                "file": f"{name}.stl",
                "triangle_count": len(artifacts[name]),
                "sha256": hashlib.sha256((output / f"{name}.stl").read_bytes()).hexdigest(),
            }
            for name in ("wall", "inlet", "outlet")
        },
    }
    # The Gmsh/CGNS conversion vocabulary uses inflow/outflow.
    region_manifest["regions"]["inflow"] = region_manifest["regions"].pop("inlet")
    region_manifest["regions"]["outflow"] = region_manifest["regions"].pop("outlet")
    region_manifest_path = output / "region_manifest.json"
    region_manifest_path.write_text(
        json.dumps(region_manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    report = {
        "operation": "generic_stl_surface_cut",
        "input_stl": str(source_path),
        "analytic_geometry_metadata_read": False,
        "cutting_planes": {
            "inlet": {
                "point": args.inlet_point, "normal": args.inlet_normal,
                "keep": args.inlet_keep,
            },
            "outlet": {
                "point": args.outlet_point, "normal": args.outlet_normal,
                "keep": args.outlet_keep,
            },
        },
        "tolerance": args.tolerance,
        "counts": {
            key: result[key] for key in (
                "source_triangle_count", "side_triangle_count",
                "inlet_loop_vertex_count", "outlet_loop_vertex_count",
                "inlet_cap_triangle_count", "outlet_cap_triangle_count",
            )
        },
        "artifacts": {**paths, "region_manifest": str(region_manifest_path)},
        "next_step": "Open closed_surface.stl in STL to Gmsh; region_manifest.json preserves wall/inflow/outflow Physical Groups.",
    }
    report_path = output / "cut_report.json"
    report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
