#!/usr/bin/env python3
"""Compile a high-level rectangle description into OneFLOW's atomic layout syntax."""

import argparse
import json
import math
from pathlib import Path

try:
    from tools.mesh_layout.layout_topology import (
        LayoutBlock,
        LayoutFace,
        LayoutLine,
        LayoutPoint,
        LayoutTopology,
    )
    from tools.mesh_layout.layout_serializer import serialize_topology
except ModuleNotFoundError:
    from layout_topology import (
        LayoutBlock,
        LayoutFace,
        LayoutLine,
        LayoutPoint,
        LayoutTopology,
    )
    from layout_serializer import serialize_topology


BOUNDARY_TYPES = {
    "none": 0,
    "interface": -1,
    "extrapolation": 1,
    "wall": 2,
    "symmetry": 3,
    "farfield": 4,
    "inflow": 5,
    "outflow": 6,
    "interior": 100,
    "overset": 1000,
}


def _boundary_value(value, line_id):
    if isinstance(value, int) and not isinstance(value, bool):
        return value
    if isinstance(value, str) and value.lower() in BOUNDARY_TYPES:
        return BOUNDARY_TYPES[value.lower()]
    raise ValueError(
        f"boundary for line {line_id} must be an integer or known boundary name"
    )


def _same_coordinate(left, right):
    return math.isclose(left, right, rel_tol=0.0, abs_tol=1.0e-12)


def _build_rectangle_topology(corners, point_counts, spacing, boundary_values):
    points = tuple(
        LayoutPoint(point_id, *coordinates)
        for point_id, coordinates in enumerate(corners, start=1)
    )
    lines = tuple(
        LayoutLine(
            line_id,
            line_id,
            line_id % 4 + 1,
            point_counts[line_id - 1],
            (spacing, -1, -1),
            boundary_values[line_id - 1],
        )
        for line_id in range(1, 5)
    )
    topology = LayoutTopology(
        points=points,
        lines=lines,
        faces=(LayoutFace(1, tuple((side, side) for side in range(1, 5))),),
        blocks=(LayoutBlock(1, ((1, 1),)),),
    )
    topology.validate()
    return topology


def _format_distribution(line):
    spacing, start_spacing, end_spacing = line.distribution
    return f"Ds({line.id})={{{spacing},{start_spacing},{end_spacing}}}"


def _serialize_topology(topology):
    return serialize_topology(
        topology,
        format_coordinate,
        _format_distribution,
    )


def compile_rectangle(spec):
    if spec.get("type") != "rectangle":
        raise ValueError("type must be 'rectangle'")

    corners = spec.get("corners")
    if not isinstance(corners, list) or len(corners) != 4:
        raise ValueError("rectangle requires exactly four corners")
    if any(not isinstance(point, list) or len(point) != 3 for point in corners):
        raise ValueError("each corner must contain three coordinates")
    for point in corners:
        if any(
            not isinstance(value, (int, float))
            or isinstance(value, bool)
            or not math.isfinite(value)
            for value in point
        ):
            raise ValueError("coordinates must be finite numbers")

    x0, y0, z0 = corners[0]
    x1, y1, z1 = corners[1]
    x2, y2, z2 = corners[2]
    x3, y3, z3 = corners[3]
    aligned = (
        _same_coordinate(x0, x3) and _same_coordinate(x1, x2)
        and _same_coordinate(y0, y1) and _same_coordinate(y2, y3)
        and _same_coordinate(z0, z1) and _same_coordinate(z0, z2)
        and _same_coordinate(z0, z3)
    )
    non_degenerate = not _same_coordinate(x0, x1) and not _same_coordinate(y0, y2)
    if not aligned or not non_degenerate:
        raise ValueError(
            "corners must describe a non-degenerate axis-aligned rectangle "
            "in perimeter order"
        )

    point_counts = spec.get("edgePointCounts")
    if not isinstance(point_counts, list) or len(point_counts) != 4:
        raise ValueError("edgePointCounts must contain four values")
    if any(not isinstance(value, int) or isinstance(value, bool) or value < 2
           for value in point_counts):
        raise ValueError("each edge point count must be an integer >= 2")
    if point_counts[0] != point_counts[2] or point_counts[1] != point_counts[3]:
        raise ValueError("opposite edges must have matching point counts for a structured rectangle")

    spacing = spec.get("spacing", "d")
    if not isinstance(spacing, str) or not spacing:
        raise ValueError("spacing must be a non-empty string")

    boundaries = spec.get("boundaries")
    if boundaries is None:
        boundary_values = [1, 2, 3, 4]
    else:
        if not isinstance(boundaries, dict):
            raise ValueError("boundaries must be an object")
        boundary_values = [
            _boundary_value(boundaries.get(f"line{line_id}", line_id), line_id)
            for line_id in range(1, 5)
        ]

    topology = _build_rectangle_topology(
        corners, point_counts, spacing, boundary_values
    )
    return _serialize_topology(topology)


def render_rectangle_svg(spec):
    corners = spec["corners"]
    point_counts = spec["edgePointCounts"]
    boundaries = spec.get("boundaries", {})
    boundary_names = {
        "none": "None", "interface": "Interface", "extrapolation": "Extrapolation",
        "wall": "Wall", "symmetry": "Symmetry", "farfield": "Farfield",
        "inflow": "Inflow", "outflow": "Outflow", "interior": "Interior",
        "overset": "Overset", 0: "None", -1: "Interface", 1: "Extrapolation",
        2: "Wall", 3: "Symmetry", 4: "Farfield", 5: "Inflow", 6: "Outflow",
        100: "Interior", 1000: "Overset",
    }

    def boundary_label(line_id):
        value = boundaries.get(f"line{line_id}", line_id)
        if isinstance(value, str):
            return boundary_names.get(value.lower(), value)
        return boundary_names.get(value, str(value))

    width, height = 1000, 640
    left_gutter, right_gutter, top_gutter, bottom_gutter = 190, 190, 105, 115
    xs = [point[0] for point in corners]
    ys = [point[1] for point in corners]
    x_span = max(max(xs) - min(xs), 1.0e-12)
    y_span = max(max(ys) - min(ys), 1.0e-12)
    scale = min(
        (width - left_gutter - right_gutter) / x_span,
        (height - top_gutter - bottom_gutter) / y_span,
    )
    draw_width, draw_height = x_span * scale, y_span * scale
    geometry_left = left_gutter + (
        width - left_gutter - right_gutter - draw_width
    ) / 2
    geometry_top = top_gutter + (
        height - top_gutter - bottom_gutter - draw_height
    ) / 2
    geometry_bottom = geometry_top + draw_height

    def screen(point):
        return (
            geometry_left + (point[0] - min(xs)) * scale,
            geometry_bottom - (point[1] - min(ys)) * scale,
        )

    screen_points = [screen(point) for point in corners]
    center_x = sum(point[0] for point in screen_points) / len(screen_points)
    center_y = sum(point[1] for point in screen_points) / len(screen_points)
    elements = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" '
        f'height="{height}" viewBox="0 0 {width} {height}">',
        '<rect x="0" y="0" width="100%" height="100%" fill="#f5f7fa"/>',
        f'<polygon points="{" ".join(f"{x:.3f},{y:.3f}" for x, y in screen_points)}" '
        'fill="#eef5ff" stroke="#202020" stroke-width="2"/>',
        f'<text x="{width / 2:.1f}" y="40" text-anchor="middle" '
        'font-size="23" fill="#202020">Rectangle layout</text>',
    ]

    for index, point in enumerate(screen_points):
        next_point = screen_points[(index + 1) % len(screen_points)]
        midpoint_x = (point[0] + next_point[0]) / 2
        midpoint_y = (point[1] + next_point[1]) / 2
        outward_x = midpoint_x - center_x
        outward_y = midpoint_y - center_y
        outward_length = (outward_x * outward_x + outward_y * outward_y) ** 0.5 or 1.0
        outward_x /= outward_length
        outward_y /= outward_length
        label_x = midpoint_x + outward_x * 34
        label_y = midpoint_y + outward_y * 34
        anchor = "middle"
        if abs(outward_x) > 0.65:
            anchor = "start" if outward_x > 0 else "end"
        elements.append(
            f'<text x="{label_x:.1f}" y="{label_y:.1f}" text-anchor="{anchor}" '
            'font-size="13" fill="#202020">'
            f'<tspan>{boundary_label(index + 1)} (Line {index + 1})</tspan>'
            f'<tspan x="{label_x:.1f}" dy="17">({point_counts[index]} points)</tspan>'
            '</text>'
        )
        elements.append(
            f'<circle cx="{point[0]:.3f}" cy="{point[1]:.3f}" r="5" fill="#d1495b"/>'
        )
        point_out_x = point[0] - center_x
        point_out_y = point[1] - center_y
        point_out_length = (point_out_x * point_out_x + point_out_y * point_out_y) ** 0.5 or 1.0
        point_label_x = point[0] + point_out_x / point_out_length * 13
        point_label_y = point[1] + point_out_y / point_out_length * 13
        elements.append(
            f'<text x="{point_label_x:.1f}" y="{point_label_y:.1f}" '
            'font-size="14" fill="#202020">'
            f'Point {index + 1}</text>'
        )

    elements.append("</svg>")
    return "\\n".join(elements) + "\\n"


def format_coordinate(value):
    if not isinstance(value, (int, float)) or isinstance(value, bool):
        raise ValueError("coordinates must be numbers")
    return format(value, ".15g")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="High-level rectangle JSON")
    parser.add_argument("-o", "--output", type=Path, help="Atomic layout output")
    parser.add_argument("--svg-output", type=Path, help="Optional topology preview SVG")
    args = parser.parse_args()

    spec = json.loads(args.input.read_text(encoding="utf-8"))
    result = compile_rectangle(spec)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(result, encoding="utf-8", newline="\n")
    else:
        print(result, end="")

    if args.svg_output:
        args.svg_output.parent.mkdir(parents=True, exist_ok=True)
        args.svg_output.write_text(
            render_rectangle_svg(spec), encoding="utf-8", newline="\n"
        )


if __name__ == "__main__":
    main()
