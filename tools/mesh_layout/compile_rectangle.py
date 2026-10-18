#!/usr/bin/env python3
"""Compile a high-level rectangle description into OneFLOW's atomic layout syntax."""

import argparse
import json
import math
from pathlib import Path


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

    lines = []
    for index, point in enumerate(corners, start=1):
        coords = ",".join(format_coordinate(value) for value in point)
        lines.append(f"Point({index})= {{{coords}}}")

    lines.append("")
    for index in range(1, 5):
        next_index = index % 4 + 1
        lines.append(f"Line({index})={{{index},{next_index}}}")

    lines.append("")
    for index, count in enumerate(point_counts, start=1):
        lines.append(f"Dim({index})={{{count}}}")

    lines.append("")
    for index in range(1, 5):
        lines.append(f"Ds({index})={{{spacing},-1,-1}}")

    lines.append("")
    for index, boundary_value in enumerate(boundary_values, start=1):
        lines.append(f"Boundary({index})={{{boundary_value}}}")

    lines.append("")
    for index in range(1, 5):
        lines.append(f"Add L2F 1 {index} {index}")

    lines.extend(["", "Add F2B 1 1 1", ""])
    return "\n".join(lines)


def render_rectangle_svg(spec):
    corners = spec["corners"]
    point_counts = spec["edgePointCounts"]
    width, height, margin = 720, 520, 80
    xs = [point[0] for point in corners]
    ys = [point[1] for point in corners]
    x_span = max(max(xs) - min(xs), 1.0e-12)
    y_span = max(max(ys) - min(ys), 1.0e-12)
    scale = min((width - 2 * margin) / x_span,
                (height - 2 * margin) / y_span)
    draw_width, draw_height = x_span * scale, y_span * scale
    left = (width - draw_width) / 2
    top = (height - draw_height) / 2

    def screen(point):
        return (left + (point[0] - min(xs)) * scale,
                height - top - (point[1] - min(ys)) * scale)

    elements = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
    ]
    for index, point in enumerate(corners):
        next_point = corners[(index + 1) % 4]
        x1, y1 = screen(point)
        x2, y2 = screen(next_point)
        mx, my = (x1 + x2) / 2, (y1 + y2) / 2
        dx, dy = x2 - x1, y2 - y1
        length = (dx * dx + dy * dy) ** 0.5 or 1.0
        nx, ny = -dy / length, dx / length
        elements.append(
            f'<line x1="{x1:.3f}" y1="{y1:.3f}" x2="{x2:.3f}" y2="{y2:.3f}" '
            f'data-line-id="{index + 1}"/>'
        )
        elements.append(
            f'<text x="{mx + nx * 24:.3f}" y="{my + ny * 24:.3f}" '
            f'text-anchor="middle" font-size="16" fill="#202020">'
            f'Line {index + 1}: {point_counts[index]} points</text>'
        )
        elements.append(
            f'<circle cx="{x1:.3f}" cy="{y1:.3f}" r="5" fill="#d1495b"/>'
        )
        elements.append(
            f'<text x="{x1 + 10:.3f}" y="{y1 - 10:.3f}" font-size="15" '
            f'fill="#202020">Point {index + 1}</text>'
        )
    elements.append(
        f'<text x="{width / 2}" y="32" text-anchor="middle" font-size="20" '
        f'fill="#202020">Rectangle topology preview</text>'
    )
    elements.append("</svg>")
    return "\n".join(elements) + "\n"


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
