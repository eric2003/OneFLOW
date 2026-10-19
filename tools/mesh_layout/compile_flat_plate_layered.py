"""Compile a two-block flat plate with explicit near-wall line segments."""

from __future__ import annotations

from tools.mesh_layout.compile_flat_plate import (
    _format_coordinate,
    _format_distribution,
    validate_flat_plate_spec,
)
from tools.mesh_layout.layout_topology import (
    LayoutBlock,
    LayoutFace,
    LayoutLine,
    LayoutPoint,
    LayoutTopology,
)
from tools.mesh_layout.layout_serializer import serialize_topology


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

    points = (
        LayoutPoint(1, x0, 0.0, 0.0),
        LayoutPoint(2, x1, 0.0, 0.0),
        LayoutPoint(3, x2, 0.0, 0.0),
        LayoutPoint(4, x0, height, 0.0),
        LayoutPoint(5, x1, height, 0.0),
        LayoutPoint(6, x2, height, 0.0),
        LayoutPoint(7, x0, y_layer, 0.0),
        LayoutPoint(8, x1, y_layer, 0.0),
        LayoutPoint(9, x2, y_layer, 0.0),
    )
    lines = (
        (1, 1, 2, base_spec["upstreamPointCount"], base_spec["distributions"]["upstream"], "upstreamBottom"),
        (2, 2, 3, base_spec["platePointCount"], base_spec["distributions"]["plate"], "plate"),
        (3, 1, 7, near_wall_points, [-1, -1], "left"),
        (4, 7, 4, upper_point_count, {"type": "tanh", "startSpacing": near_wall_distance / near_wall_intervals, "endSpacing": "auto"}, "left"),
        (5, 2, 8, near_wall_points, [-1, -1], "interface"),
        (6, 8, 5, upper_point_count, {"type": "tanh", "startSpacing": near_wall_distance / near_wall_intervals, "endSpacing": "auto"}, "interface"),
        (7, 3, 9, near_wall_points, [-1, -1], "right"),
        (8, 9, 6, upper_point_count, {"type": "tanh", "startSpacing": near_wall_distance / near_wall_intervals, "endSpacing": "auto"}, "right"),
        (9, 4, 5, base_spec["upstreamPointCount"], base_spec["distributions"].get("topUpstream", base_spec["distributions"].get("top")), "topUpstream"),
        (10, 5, 6, base_spec["platePointCount"], base_spec["distributions"].get("topPlate", base_spec["distributions"].get("top")), "topPlate"),
    )

    boundaries = base_spec["boundaries"]
    from tools.mesh_layout.compile_flat_plate import _boundary_type

    topology = LayoutTopology(
        points=points,
        lines=tuple(
            LayoutLine(
                line_id,
                start,
                end,
                count,
                distribution,
                _boundary_type(boundaries, boundary_name),
            )
            for line_id, start, end, count, distribution, boundary_name in lines
        ),
        faces=(
            LayoutFace(1, ((1, 1), (2, 5), (2, 6), (3, 9), (4, 3), (4, 4))),
            LayoutFace(2, ((1, 2), (2, 7), (2, 8), (3, 10), (4, 5), (4, 6))),
        ),
        blocks=(
            LayoutBlock(1, ((1, 1),)),
            LayoutBlock(2, ((1, 2),)),
        ),
    )
    topology.validate()
    return serialize_topology(topology, _format_coordinate, _format_distribution)


__all__ = ["compile_flat_plate_layered"]
