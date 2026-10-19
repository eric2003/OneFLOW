"""Compile a two-block flat plate with explicit near-wall line segments."""

from __future__ import annotations

from types import SimpleNamespace

from tools.mesh_layout.compile_flat_plate import (
    _boundary_type,
    _format_coordinate,
    _format_distribution,
    validate_flat_plate_spec,
)


def compile_flat_plate_layered(spec):
    """Compile a ten-line topology with 41 points in the first 0.02 near the wall."""
    base_spec = dict(spec)
    base_spec["type"] = "flat_plate"
    validate_flat_plate_spec(base_spec)

    vertical = base_spec["distributions"]["vertical"]
    if vertical.get("type") != "boundary_layer":
        raise ValueError("vertical distribution must use boundary_layer")
    near_wall_intervals = vertical["nearWallIntervalCount"]
    near_wall_distance = float(vertical["nearWallDistance"])
    normal_count = base_spec["normalPointCount"]
    near_wall_points = near_wall_intervals + 1
    if normal_count <= near_wall_intervals:
        raise ValueError("normalPointCount must exceed nearWallIntervalCount")
    upper_point_count = normal_count - near_wall_intervals
    if near_wall_distance >= float(base_spec["domainHeight"]):
        raise ValueError("nearWallDistance must be smaller than domainHeight")

    length = float(base_spec["plateLength"])
    upstream = float(base_spec["upstreamDistance"])
    height = float(base_spec["domainHeight"])
    x0, x1, x2 = -upstream, 0.0, length
    y_layer = near_wall_distance
    points = [
        (1, x0, 0.0, 0.0),
        (2, x1, 0.0, 0.0),
        (3, x2, 0.0, 0.0),
        (4, x0, height, 0.0),
        (5, x1, height, 0.0),
        (6, x2, height, 0.0),
        (7, x0, y_layer, 0.0),
        (8, x1, y_layer, 0.0),
        (9, x2, y_layer, 0.0),
    ]
    upper_distribution = {
        "type": "tanh",
        "startSpacing": near_wall_distance / near_wall_intervals,
        "endSpacing": "auto",
    }
    distributions = base_spec["distributions"]
    top_upstream = distributions.get("topUpstream", distributions.get("top"))
    top_plate = distributions.get("topPlate")
    if top_plate is None:
        if isinstance(top_upstream, list):
            top_plate = list(reversed(top_upstream))
        elif isinstance(top_upstream, dict) and top_upstream.get("type") == "tanh":
            top_plate = {
                "type": "tanh",
                "startSpacing": top_upstream["endSpacing"],
                "endSpacing": top_upstream["startSpacing"],
            }
        else:
            top_plate = top_upstream

    line_definitions = [
        (1, 1, 2, base_spec["upstreamPointCount"], distributions["upstream"], "upstreamBottom"),
        (2, 2, 3, base_spec["platePointCount"], distributions["plate"], "plate"),
        (3, 1, 7, near_wall_points, [-1, -1], "left"),
        (4, 7, 4, upper_point_count, upper_distribution, "left"),
        (5, 2, 8, near_wall_points, [-1, -1], "interface"),
        (6, 8, 5, upper_point_count, upper_distribution, "interface"),
        (7, 3, 9, near_wall_points, [-1, -1], "right"),
        (8, 9, 6, upper_point_count, upper_distribution, "right"),
        (9, 4, 5, base_spec["upstreamPointCount"], top_upstream, "topUpstream"),
        (10, 5, 6, base_spec["platePointCount"], top_plate, "topPlate"),
    ]
    output = []
    for point_id, x, y, z in points:
        coordinates = ",".join(_format_coordinate(value) for value in (x, y, z))
        output.append(f"Point({point_id})= {{{coordinates}}}")
    output.append("")
    for line_id, start, end, _, _, _ in line_definitions:
        output.append(f"Line({line_id})={{{start},{end}}}")
    output.append("")
    for line_id, _, _, count, _, _ in line_definitions:
        output.append(f"Dim({line_id})={{{count}}}")
    output.append("")
    for line_id, _, _, _, distribution, _ in line_definitions:
        output.append(_format_distribution(SimpleNamespace(id=line_id, distribution=distribution)))
    output.append("")
    for line_id, _, _, _, _, boundary_name in line_definitions:
        output.append(f"Boundary({line_id})={{{_boundary_type(base_spec['boundaries'], boundary_name)}}}")
    output.extend([
        "",
        "Add L2F 1 1 1",
        "Add L2F 1 2 5",
        "Add L2F 1 2 6",
        "Add L2F 1 3 9",
        "Add L2F 1 4 3",
        "Add L2F 1 4 4",
        "",
        "Add L2F 2 1 2",
        "Add L2F 2 2 7",
        "Add L2F 2 2 8",
        "Add L2F 2 3 10",
        "Add L2F 2 4 5",
        "Add L2F 2 4 6",
        "",
        "Add F2B 1 1 1",
        "Add F2B 2 1 2",
        "",
    ])
    return "\n".join(output)

def render_flat_plate_layered_svg(spec):
    from tools.mesh_layout.compile_flat_plate import render_flat_plate_svg

    preview_spec = dict(spec)
    preview_spec["type"] = "flat_plate"
    svg = render_flat_plate_svg(preview_spec)
    length = float(spec["plateLength"])
    upstream = float(spec["upstreamDistance"])
    height = float(spec["domainHeight"])
    vertical = spec["distributions"]["vertical"]
    distance = float(vertical["nearWallDistance"])
    near_wall_points = int(vertical["nearWallIntervalCount"]) + 1
    width, canvas_height = 1000, 640
    left_gutter, right_gutter, top_gutter, bottom_gutter = 210, 170, 95, 110
    x_min, x_max = -upstream, length
    x_span = max(x_max - x_min, 1.0e-12)
    y_span = max(height, 1.0e-12)
    scale = min(
        (width - left_gutter - right_gutter) / x_span,
        (canvas_height - top_gutter - bottom_gutter) / y_span,
    )
    draw_width = x_span * scale
    draw_height = y_span * scale
    geometry_left = left_gutter + (
        width - left_gutter - right_gutter - draw_width
    ) / 2
    geometry_top = top_gutter + (
        canvas_height - top_gutter - bottom_gutter - draw_height
    ) / 2
    layer_y = geometry_top + draw_height - distance * scale
    layer = (
        f'<line x1="{geometry_left:.3f}" y1="{layer_y:.3f}" '
        f'x2="{geometry_left + draw_width:.3f}" y2="{layer_y:.3f}" '
        'stroke="#1769aa" stroke-width="2" stroke-dasharray="5,3"/>'
        f'<text x="{geometry_left + draw_width + 8:.3f}" y="{layer_y + 4:.3f}" '
        'font-size="12" fill="#1769aa">y = near-wall distance ({near_wall_points} points)</text>'
    )
    return svg.replace("</svg>", layer + "</svg>")


__all__ = ["compile_flat_plate_layered", "render_flat_plate_layered_svg"]
