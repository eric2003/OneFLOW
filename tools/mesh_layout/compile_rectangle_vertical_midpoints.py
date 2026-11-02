"""Compile a rectangle with each vertical side split at its midpoint."""

from __future__ import annotations

import math


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


def _boundary_value(value, side):
    if isinstance(value, int) and not isinstance(value, bool):
        return value
    if isinstance(value, str) and value.lower() in BOUNDARY_TYPES:
        return BOUNDARY_TYPES[value.lower()]
    raise ValueError(f"boundary for {side} must be an integer or known boundary name")


def _number(value):
    if not isinstance(value, (int, float)) or isinstance(value, bool):
        raise ValueError("coordinates must be numbers")
    if not math.isfinite(value):
        raise ValueError("coordinates must be finite")
    return format(value, ".15g")


def compile_rectangle_vertical_midpoints(spec):
    if spec.get("type") != "rectangle_vertical_midpoints":
        raise ValueError("type must be 'rectangle_vertical_midpoints'")

    corners = spec.get("corners")
    if not isinstance(corners, list) or len(corners) != 4:
        raise ValueError("rectangle requires exactly four corners")
    if any(not isinstance(point, list) or len(point) != 3 for point in corners):
        raise ValueError("each corner must contain three coordinates")
    for point in corners:
        for value in point:
            _number(value)

    (x0, y0, z0), (x1, y1, z1), (x2, y2, z2), (x3, y3, z3) = corners
    close = lambda a, b: math.isclose(a, b, rel_tol=0.0, abs_tol=1.0e-12)
    aligned = (
        close(x0, x3) and close(x1, x2)
        and close(y0, y1) and close(y2, y3)
        and close(z0, z1) and close(z0, z2) and close(z0, z3)
    )
    if not aligned or close(x0, x1) or close(y0, y2):
        raise ValueError("corners must describe a non-degenerate axis-aligned rectangle")

    counts = spec.get("segmentPointCounts")
    if not isinstance(counts, list) or len(counts) != 6:
        raise ValueError("segmentPointCounts must contain six values")
    if any(not isinstance(value, int) or isinstance(value, bool) or value < 2 for value in counts):
        raise ValueError("each segment point count must be an integer >= 2")
    if counts[0] != counts[3] or len(set(counts[1:3] + counts[4:6])) != 1:
        raise ValueError("opposite sides and the four vertical segments must have matching point counts")

    spacing = spec.get("spacing", "d")
    if not isinstance(spacing, str) or not spacing:
        raise ValueError("spacing must be a non-empty string")

    boundaries = spec.get("boundaries", {})
    if not isinstance(boundaries, dict):
        raise ValueError("boundaries must be an object")
    boundary_values = [
        _boundary_value(boundaries.get("bottom", "extrapolation"), "bottom"),
        _boundary_value(boundaries.get("right", "wall"), "right"),
        _boundary_value(boundaries.get("right", "wall"), "right"),
        _boundary_value(boundaries.get("top", "symmetry"), "top"),
        _boundary_value(boundaries.get("left", "farfield"), "left"),
        _boundary_value(boundaries.get("left", "farfield"), "left"),
    ]

    points = [
        (x0, y0, z0), (x1, y1, z1), (x2, y2, z2), (x3, y3, z3),
        ((x1 + x2) / 2, (y1 + y2) / 2, (z1 + z2) / 2),
        ((x0 + x3) / 2, (y0 + y3) / 2, (z0 + z3) / 2),
    ]
    endpoints = [(1, 2), (2, 5), (5, 3), (3, 4), (4, 6), (6, 1)]
    lines = [
        f"Line({line_id}) = {{{start},{end}}}"
        for line_id, (start, end) in enumerate(endpoints, start=1)
    ]
    output = [f"Point({i}) = {{{','.join(_number(v) for v in point)}}}" for i, point in enumerate(points, start=1)]
    output += ["", *lines, ""]
    output += [f"Dim({i}) = {{{count}}}" for i, count in enumerate(counts, start=1)]
    output += ["", *[f"Ds({i}) = {{{spacing},-1,-1}}" for i in range(1, 7)], ""]
    output += [f"Boundary({i}) = {{{value}}}" for i, value in enumerate(boundary_values, start=1)]
    output += [
        "",
        "Add L2F 1 1 1",
        "Add L2F 1 2 2",
        "Add L2F 1 2 3",
        "Add L2F 1 3 4",
        "Add L2F 1 4 5",
        "Add L2F 1 4 6",
        "",
        "Add F2B 1 1 1",
    ]
    return "\n".join(output) + "\n"


def render_rectangle_vertical_midpoints_svg(spec):
    corners = spec["corners"]
    points = [
        corners[0], corners[1], corners[2], corners[3],
        [(corners[1][axis] + corners[2][axis]) / 2 for axis in range(3)],
        [(corners[0][axis] + corners[3][axis]) / 2 for axis in range(3)],
    ]
    width, height = 900, 620
    margin = 150
    xs = [point[0] for point in corners]
    ys = [point[1] for point in corners]
    x_span = max(max(xs) - min(xs), 1.0e-12)
    y_span = max(max(ys) - min(ys), 1.0e-12)
    scale = min((width - 2 * margin) / x_span, (height - 2 * margin) / y_span)

    def screen(point):
        return (margin + (point[0] - min(xs)) * scale,
                height - margin - (point[1] - min(ys)) * scale)

    screen_points = [screen(point) for point in points]
    line_pairs = [(0, 1), (1, 4), (4, 2), (2, 3), (3, 5), (5, 0)]
    boundaries = spec.get("boundaries", {})
    labels = [
        str(boundaries.get("bottom", "extrapolation")),
        str(boundaries.get("right", "wall")),
        str(boundaries.get("right", "wall")),
        str(boundaries.get("top", "symmetry")),
        str(boundaries.get("left", "farfield")),
        str(boundaries.get("left", "farfield")),
    ]
    elements = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#f5f7fa"/>',
        '<text x="450" y="40" text-anchor="middle" font-size="22" fill="#202020">Rectangle with vertical midpoints</text>',
    ]
    for index, (start, end) in enumerate(line_pairs, start=1):
        x1, y1 = screen_points[start]
        x2, y2 = screen_points[end]
        elements.append(f'<line x1="{x1:.3f}" y1="{y1:.3f}" x2="{x2:.3f}" y2="{y2:.3f}" stroke="#202020" stroke-width="2"/>')
        elements.append(f'<text x="{(x1+x2)/2:.3f}" y="{(y1+y2)/2-10:.3f}" font-size="12" text-anchor="middle" fill="#202020">Line {index}: {labels[index-1]}</text>')
    for index, (x, y) in enumerate(screen_points, start=1):
        elements.append(f'<circle cx="{x:.3f}" cy="{y:.3f}" r="5" fill="#d1495b"/>')
        elements.append(f'<text x="{x+8:.3f}" y="{y-8:.3f}" font-size="13" fill="#202020">Point {index}</text>')
    elements.append('</svg>')
    return "\n".join(elements) + "\n"


__all__ = ["compile_rectangle_vertical_midpoints", "render_rectangle_vertical_midpoints_svg"]
