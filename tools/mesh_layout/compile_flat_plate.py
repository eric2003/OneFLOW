#!/usr/bin/env python3
"""Compile a flat-plate boundary-layer specification into OneFLOW layout syntax."""

import argparse
import json
from pathlib import Path

try:
    from tools.mesh_layout.layout_topology import (
        LayoutBlock,
        LayoutFace,
        LayoutLine,
        LayoutPoint,
        LayoutTopology,
    )
except ModuleNotFoundError:
    from layout_topology import (
        LayoutBlock,
        LayoutFace,
        LayoutLine,
        LayoutPoint,
        LayoutTopology,
    )


BC_TYPES = {
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


def _require_positive_int(spec, key, minimum=2):
    value = spec.get(key)
    if not isinstance(value, int) or isinstance(value, bool) or value < minimum:
        raise ValueError(f"{key} must be an integer >= {minimum}")
    return value


def _require_number(spec, key):
    value = spec.get(key)
    if not isinstance(value, (int, float)) or isinstance(value, bool):
        raise ValueError(f"{key} must be a number")
    return float(value)


def _boundary_type(boundaries, name):
    value = boundaries.get(name)
    if isinstance(value, int) and not isinstance(value, bool):
        return value
    if isinstance(value, str) and value.lower() in BC_TYPES:
        return BC_TYPES[value.lower()]
    raise ValueError(f"boundary {name} must be an integer or known boundary name")


def _format_coordinate(value):
    return format(value, ".15g")


def _format_spacing(value):
    if value == "auto":
        return "-1"
    formatted = format(float(value), ".15g")
    if "e" in formatted.lower():
        mantissa, exponent = formatted.lower().split("e", 1)
        if "." not in mantissa:
            mantissa += ".0"
        formatted = f"{mantissa}e{int(exponent)}"
    return formatted



def validate_flat_plate_spec(spec):
    if spec.get("type") != "flat_plate":
        raise ValueError("type must be 'flat_plate'")

    length = _require_number(spec, "plateLength")
    upstream = _require_number(spec, "upstreamDistance")
    height = _require_number(spec, "domainHeight")
    if length <= 0 or upstream <= 0 or height <= 0:
        raise ValueError("plateLength, upstreamDistance, and domainHeight must be positive")

    _require_positive_int(spec, "upstreamPointCount")
    _require_positive_int(spec, "platePointCount")
    _require_positive_int(spec, "normalPointCount")

    distributions = spec.get("distributions", {})
    if not isinstance(distributions, dict):
        raise ValueError("distributions must be an object")
    names = ["upstream", "plate", "vertical"]
    has_top_upstream = "topUpstream" in distributions
    has_top_plate = "topPlate" in distributions
    if has_top_upstream != has_top_plate:
        raise ValueError("topUpstream and topPlate distributions must be specified together")
    names.extend(["topUpstream", "topPlate"] if has_top_upstream else ["top"])
    for name in names:
        value = distributions.get(name)
        if isinstance(value, list):
            if len(value) != 2:
                raise ValueError(f"distribution {name} must contain two values")
        elif isinstance(value, dict):
            distribution_type = value.get("type")
            if distribution_type == "boundary_layer":
                count = value.get("nearWallPointCount")
                distance = value.get("nearWallDistance")
                if not isinstance(count, int) or isinstance(count, bool) or count < 2:
                    raise ValueError(f"distribution {name} nearWallPointCount must be an integer >= 2")
                if not isinstance(distance, (int, float)) or isinstance(distance, bool) or distance <= 0:
                    raise ValueError(f"distribution {name} nearWallDistance must be positive")
            elif distribution_type == "tanh":
                for endpoint in ("startSpacing", "endSpacing"):
                    spacing = value.get(endpoint)
                    if spacing == "auto":
                        continue
                    if not isinstance(spacing, (int, float)) or isinstance(spacing, bool) or spacing <= 0:
                        raise ValueError(f"distribution {name} {endpoint} must be positive or auto")
            else:
                raise ValueError(f"distribution {name} has unsupported object type")
        else:
            raise ValueError(f"distribution {name} must contain two values or a supported distribution object")

    boundaries = spec.get("boundaries")
    if not isinstance(boundaries, dict):
        raise ValueError("boundaries must be an object")
    for name in ("upstreamBottom", "plate", "left", "interface", "right",
                 "topUpstream", "topPlate"):
        _boundary_type(boundaries, name)

    return spec


def compile_flat_plate(spec):
    validate_flat_plate_spec(spec)

    length = float(spec["plateLength"])
    upstream = float(spec["upstreamDistance"])
    height = float(spec["domainHeight"])
    upstream_count = spec["upstreamPointCount"]
    plate_count = spec["platePointCount"]
    normal_count = spec["normalPointCount"]
    x0 = -upstream
    x1 = 0.0
    x2 = length

    points = [
        (1, x0, 0.0, 0.0),
        (2, x1, 0.0, 0.0),
        (3, x2, 0.0, 0.0),
        (4, x0, height, 0.0),
        (5, x1, height, 0.0),
        (6, x2, height, 0.0),
    ]

    lines = [
        (1, 1, 2),
        (2, 2, 3),
        (3, 1, 4),
        (4, 2, 5),
        (5, 3, 6),
        (6, 4, 5),
        (7, 5, 6),
    ]

    counts = [
        upstream_count,
        plate_count,
        normal_count,
        normal_count,
        normal_count,
        upstream_count,
        plate_count,
    ]

    distributions = spec["distributions"]
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

    distribution_values = [
        distributions["upstream"],
        distributions["plate"],
        distributions["vertical"],
        distributions["vertical"],
        distributions["vertical"],
        top_upstream,
        top_plate,
    ]

    boundary_values = [
        _boundary_type(spec["boundaries"], "upstreamBottom"),
        _boundary_type(spec["boundaries"], "plate"),
        _boundary_type(spec["boundaries"], "left"),
        _boundary_type(spec["boundaries"], "interface"),
        _boundary_type(spec["boundaries"], "right"),
        _boundary_type(spec["boundaries"], "topUpstream"),
        _boundary_type(spec["boundaries"], "topPlate"),
    ]

    topology = LayoutTopology(
        points=tuple(LayoutPoint(point_id, x, y, z) for point_id, x, y, z in points),
        lines=tuple(
            LayoutLine(
                line_id, p1, p2, counts[line_id - 1],
                distribution_values[line_id - 1], boundary_values[line_id - 1],
            )
            for line_id, p1, p2 in lines
        ),
        faces=(
            LayoutFace(1, ((1, 1), (2, 4), (3, 6), (4, 3))),
            LayoutFace(2, ((1, 2), (2, 5), (3, 7), (4, 4))),
        ),
        blocks=(
            LayoutBlock(1, ((1, 1),)),
            LayoutBlock(2, ((1, 2),)),
        ),
    )
    topology.validate()
    lines_out = []
    for point in topology.points:
        lines_out.append(
            f"Point({point.id})= {{{_format_coordinate(point.x)},"
            f"{_format_coordinate(point.y)},{_format_coordinate(point.z)}}}"
        )

    lines_out.append("")
    for line in topology.lines:
        lines_out.append(
            f"Line({line.id})={{{line.start_point},{line.end_point}}}"
        )

    lines_out.append("")
    for line in topology.lines:
        lines_out.append(f"Dim({line.id})={{{line.point_count}}}")

    lines_out.append("")
    for line in topology.lines:
        line_id = line.id
        values = line.distribution
        if isinstance(values, list):
            start_value, end_value = values
            lines_out.append(f"Ds({line_id})={{d,{start_value},{end_value}}}")
        elif values["type"] == "boundary_layer":
            lines_out.append(
                f"Ds({line_id})={{bl,{values['nearWallPointCount']},"
                f"{values['nearWallDistance']}}}"
            )
        else:
            start_value = _format_spacing(values["startSpacing"])
            end_value = _format_spacing(values["endSpacing"])
            lines_out.append(f"Ds({line_id})={{tanh,{start_value},{end_value}}}")

    lines_out.append("")
    for line in topology.lines:
        lines_out.append(f"Boundary({line.id})={{{line.boundary_type}}}")

    for face in topology.faces:
        lines_out.append("")
        for position, line_id in face.line_by_side:
            lines_out.append(f"Add L2F {face.id} {position} {line_id}")

    lines_out.append("")
    for block in topology.blocks:
        for position, face_id in block.face_by_position:
            lines_out.append(f"Add F2B {block.id} {position} {face_id}")
    lines_out.append("")
    return "\n".join(lines_out)


def render_flat_plate_svg(spec):
    validate_flat_plate_spec(spec)

    length = float(spec["plateLength"])
    upstream = float(spec["upstreamDistance"])
    height = float(spec["domainHeight"])
    x_min, x_max = -upstream, length

    width, canvas_height, margin = 900, 560, 80
    x_span = max(x_max - x_min, 1.0e-12)
    y_span = max(height, 1.0e-12)
    scale = min((width - 2 * margin) / x_span,
                (canvas_height - 2 * margin) / y_span)

    draw_width = x_span * scale
    draw_height = y_span * scale
    left = (width - draw_width) / 2
    bottom = canvas_height - (canvas_height - draw_height) / 2

    def screen(x, y):
        return left + (x - x_min) * scale, bottom - y * scale

    p1 = screen(-upstream, 0.0)
    p2 = screen(0.0, 0.0)
    p3 = screen(length, 0.0)
    p4 = screen(-upstream, height)
    p5 = screen(0.0, height)
    p6 = screen(length, height)

    boundaries = spec["boundaries"]
    boundary_names = {
        "none": "None", "interface": "Interface", "extrapolation": "Extrapolation",
        "wall": "Wall", "symmetry": "Symmetry", "farfield": "Farfield",
        "inflow": "Inflow", "outflow": "Outflow", "interior": "Interior",
        "overset": "Overset", 0: "None", -1: "Interface", 1: "Extrapolation",
        2: "Wall", 3: "Symmetry", 4: "Farfield", 5: "Inflow", 6: "Outflow",
        100: "Interior", 1000: "Overset",
    }

    def boundary_label(name):
        value = boundaries[name]
        return boundary_names.get(value.lower(), value) if isinstance(value, str) else boundary_names.get(value, str(value))

    edge_labels = [
        ("upstreamBottom", boundary_label("upstreamBottom"), spec["upstreamPointCount"], (p1, p2), "bottom"),
        ("plate", boundary_label("plate"), spec["platePointCount"], (p2, p3), "bottom"),
        ("left", boundary_label("left"), spec["normalPointCount"], (p1, p4), "left"),
        ("interface", boundary_label("interface"), spec["normalPointCount"], (p2, p5), "right"),
        ("right", boundary_label("right"), spec["normalPointCount"], (p3, p6), "right"),
        ("topUpstream", boundary_label("topUpstream"), spec["upstreamPointCount"], (p4, p5), "top"),
        ("topPlate", boundary_label("topPlate"), spec["platePointCount"], (p5, p6), "top"),
    ]

    elements = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" '
        f'height="{canvas_height}" viewBox="0 0 {width} {canvas_height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        f'<rect x="{p1[0]:.3f}" y="{p4[1]:.3f}" width="{p6[0]-p1[0]:.3f}" '
        f'height="{p1[1]-p4[1]:.3f}" fill="#eef5ff" stroke="#202020"/>',
        f'<line x1="{p2[0]:.3f}" y1="{p2[1]:.3f}" x2="{p5[0]:.3f}" '
        f'y2="{p5[1]:.3f}" stroke="#888" stroke-dasharray="6,5"/>',
        f'<line x1="{p2[0]:.3f}" y1="{p2[1]:.3f}" x2="{p3[0]:.3f}" '
        f'y2="{p3[1]:.3f}" stroke="#b23a48" stroke-width="5"/>',
        f'<text x="{width/2:.1f}" y="34" text-anchor="middle" '
        f'font-size="21" fill="#202020">Flat plate boundary-layer layout</text>',
    ]
    for edge_name, boundary_name, point_count, endpoints, side in edge_labels:
        (x1, y1), (x2, y2) = endpoints
        if side == "bottom":
            x = (x1 + x2) / 2
            y = y1 + (55 if edge_name == "plate" else 34)
            anchor = "middle"
        elif side == "top":
            x = (x1 + x2) / 2
            y = y1 - 14
            anchor = "middle"
        elif side == "left":
            x = x1 - 8
            y = (y1 + y2) / 2
            anchor = "end"
        else:
            x = x1 + 8
            y = (y1 + y2) / 2
            anchor = "start"
        elements.append(
            f'<text x="{x:.1f}" y="{y:.1f}" text-anchor="{anchor}" '
            f'font-size="12" fill="#202020">{edge_name}: {boundary_name} '
            f'({point_count} points)</text>'
        )
    elements.append("</svg>")
    return "\n".join(elements) + "\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="Flat-plate JSON specification")
    parser.add_argument("-o", "--output", type=Path, help="Atomic layout output")
    parser.add_argument("--svg-output", type=Path, help="Optional topology preview SVG")
    args = parser.parse_args()

    spec = json.loads(args.input.read_text(encoding="utf-8"))
    output = compile_flat_plate(spec)

    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(output, encoding="utf-8", newline="\n")
    else:
        print(output, end="")

    if args.svg_output:
        args.svg_output.parent.mkdir(parents=True, exist_ok=True)
        args.svg_output.write_text(
            render_flat_plate_svg(spec), encoding="utf-8", newline="\n"
        )


if __name__ == "__main__":
    main()
