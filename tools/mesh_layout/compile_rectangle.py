#!/usr/bin/env python3
"""Compile a high-level rectangle description into OneFLOW's atomic layout syntax."""

import argparse
import json
from pathlib import Path


def compile_rectangle(spec):
    if spec.get("type") != "rectangle":
        raise ValueError("type must be 'rectangle'")

    corners = spec.get("corners")
    if not isinstance(corners, list) or len(corners) != 4:
        raise ValueError("rectangle requires exactly four corners")
    if any(not isinstance(point, list) or len(point) != 3 for point in corners):
        raise ValueError("each corner must contain three coordinates")

    point_counts = spec.get("edgePointCounts")
    if not isinstance(point_counts, list) or len(point_counts) != 4:
        raise ValueError("edgePointCounts must contain four values")
    if any(not isinstance(value, int) or isinstance(value, bool) or value < 2
           for value in point_counts):
        raise ValueError("each edge point count must be an integer >= 2")

    spacing = spec.get("spacing", "d")
    if not isinstance(spacing, str) or not spacing:
        raise ValueError("spacing must be a non-empty string")

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
    for index in range(1, 5):
        lines.append(f"Boundary({index})={{{index}}}")

    lines.append("")
    for index in range(1, 5):
        lines.append(f"Add L2F 1 {index} {index}")

    lines.extend(["", "Add F2B 1 1 1", ""])
    return "\n".join(lines)


def format_coordinate(value):
    if not isinstance(value, (int, float)) or isinstance(value, bool):
        raise ValueError("coordinates must be numbers")
    return format(value, ".15g")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="High-level rectangle JSON")
    parser.add_argument("-o", "--output", type=Path, help="Atomic layout output")
    args = parser.parse_args()

    spec = json.loads(args.input.read_text(encoding="utf-8"))
    result = compile_rectangle(spec)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(result, encoding="utf-8", newline="\n")
    else:
        print(result, end="")


if __name__ == "__main__":
    main()
