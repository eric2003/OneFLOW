#!/usr/bin/env python3
"""Export a single-zone OneFLOW structured-grid file for inspection.

The current case reads the binary header and coordinate arrays written by
GridMachine: zone count, ni/nj/nk, then x/y/z arrays for each zone.
The exporter intentionally supports one 2D zone first; it rejects other
layouts rather than guessing their topology.
"""

import argparse
import json
import math
import struct
from pathlib import Path


def read_int(data, offset, endian):
    size = 4
    if offset + size > len(data):
        raise ValueError("Unexpected end of file while reading a 32-bit integer")
    return struct.unpack_from(endian + "i", data, offset)[0], offset + size


def read_real_array(data, offset, count, endian, real_bytes):
    code = "d" if real_bytes == 8 else "f"
    size = struct.calcsize(code) * count
    if offset + size > len(data):
        raise ValueError("Unexpected end of file while reading coordinates")
    values = list(struct.unpack_from(endian + str(count) + code, data, offset))
    if not all(math.isfinite(value) for value in values):
        raise ValueError("Coordinate array contains a non-finite value")
    return values, offset + size


def read_structured_grid(path, endian, real_bytes):
    data = path.read_bytes()
    offset = 0
    zone_count, offset = read_int(data, offset, endian)
    if zone_count != 1:
        raise ValueError(
            f"This first exporter supports one zone; input contains {zone_count}"
        )

    ni, offset = read_int(data, offset, endian)
    nj, offset = read_int(data, offset, endian)
    nk, offset = read_int(data, offset, endian)
    if ni < 2 or nj < 2 or nk != 1:
        raise ValueError(
            f"Expected a 2D structured grid with ni,nj >= 2 and nk == 1; "
            f"got ({ni}, {nj}, {nk})"
        )

    node_count = ni * nj * nk
    x, offset = read_real_array(data, offset, node_count, endian, real_bytes)
    y, offset = read_real_array(data, offset, node_count, endian, real_bytes)
    z, offset = read_real_array(data, offset, node_count, endian, real_bytes)
    points = [[x[i], y[i], z[i]] for i in range(node_count)]
    return ni, nj, points


def build_topology(ni, nj):
    cells = []
    face_map = {}
    faces = []

    def add_face(node_a, node_b, cell_id):
        key = tuple(sorted((node_a, node_b)))
        if key in face_map:
            faces[face_map[key]]["neighbour"] = cell_id
        else:
            face_id = len(faces)
            face_map[key] = face_id
            faces.append({
                "id": face_id,
                "nodes": [node_a, node_b],
                "owner": cell_id,
                "neighbour": -1,
            })

    for j in range(nj - 1):
        for i in range(ni - 1):
            n0 = j * ni + i
            n1 = n0 + 1
            n3 = (j + 1) * ni + i
            n2 = n3 + 1
            cell_id = len(cells)
            cells.append({"id": cell_id, "type": "quad", "nodes": [n0, n1, n2, n3]})
            add_face(n0, n1, cell_id)
            add_face(n1, n2, cell_id)
            add_face(n2, n3, cell_id)
            add_face(n3, n0, cell_id)

    boundary_faces = [
        {"face": face["id"], "region": "boundary"}
        for face in faces if face["neighbour"] == -1
    ]
    return cells, faces, boundary_faces


def coordinate_bounds(points):
    if not points:
        raise ValueError("Cannot compute coordinate bounds for an empty mesh")
    return {
        "min": [min(point[axis] for point in points) for axis in range(3)],
        "max": [max(point[axis] for point in points) for axis in range(3)],
    }


def write_json(path, document):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(document, indent=2) + "\n", encoding="utf-8")


def write_vtu(path, points, cells):
    path.parent.mkdir(parents=True, exist_ok=True)
    connectivity = " ".join(str(node) for cell in cells for node in cell["nodes"])
    offsets = " ".join(str(4 * (i + 1)) for i in range(len(cells)))
    types = " ".join("9" for _ in cells)  # VTK_QUAD
    point_values = "\n".join(" ".join(f"{v:.17g}" for v in point) for point in points)
    content = f"""<?xml version="1.0"?>
<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">
  <UnstructuredGrid>
    <Piece NumberOfPoints="{len(points)}" NumberOfCells="{len(cells)}">
      <Points>
        <DataArray type="Float64" NumberOfComponents="3" format="ascii">
{point_values}
        </DataArray>
      </Points>
      <Cells>
        <DataArray type="Int64" Name="connectivity" format="ascii">{connectivity}</DataArray>
        <DataArray type="Int64" Name="offsets" format="ascii">{offsets}</DataArray>
        <DataArray type="UInt8" Name="types" format="ascii">{types}</DataArray>
      </Cells>
    </Piece>
  </UnstructuredGrid>
</VTKFile>
"""
    path.write_text(content, encoding="utf-8")


def write_svg(path, ni, nj, points):
    path.parent.mkdir(parents=True, exist_ok=True)
    xs = [point[0] for point in points]
    ys = [point[1] for point in points]
    xmin, xmax = min(xs), max(xs)
    ymin, ymax = min(ys), max(ys)
    width, height, margin = 800, 600, 40
    dx = xmax - xmin or 1.0
    dy = ymax - ymin or 1.0
    scale = min((width - 2 * margin) / dx, (height - 2 * margin) / dy)
    draw_width, draw_height = dx * scale, dy * scale
    left = (width - draw_width) / 2
    top = (height - draw_height) / 2

    def screen(point):
        return (left + (point[0] - xmin) * scale,
                height - top - (point[1] - ymin) * scale)

    lines = []
    for j in range(nj):
        for i in range(ni - 1):
            a = screen(points[j * ni + i])
            b = screen(points[j * ni + i + 1])
            lines.append(f'<line x1="{a[0]:.3f}" y1="{a[1]:.3f}" '
                         f'x2="{b[0]:.3f}" y2="{b[1]:.3f}"/>')
    for i in range(ni):
        for j in range(nj - 1):
            a = screen(points[j * ni + i])
            b = screen(points[(j + 1) * ni + i])
            lines.append(f'<line x1="{a[0]:.3f}" y1="{a[1]:.3f}" '
                         f'x2="{b[0]:.3f}" y2="{b[1]:.3f}"/>')
    markers = []
    for node_id, point in enumerate(points):
        x, y = screen(point)
        markers.append(
            f'<circle cx="{x:.3f}" cy="{y:.3f}" r="3" fill="#d1495b" '
            f'data-node-id="{node_id}"/>'
        )

    cell_labels = []
    for j in range(nj - 1):
        for i in range(ni - 1):
            node_ids = (j * ni + i, j * ni + i + 1,
                        (j + 1) * ni + i + 1, (j + 1) * ni + i)
            center = [sum(points[node_id][axis] for node_id in node_ids) / 4
                      for axis in range(3)]
            x, y = screen(center)
            cell_id = j * (ni - 1) + i
            cell_labels.append(
                f'<text x="{x:.3f}" y="{y:.3f}" text-anchor="middle" '
                f'dominant-baseline="middle" font-size="12" fill="#202020" '
                f'data-cell-id="{cell_id}">{cell_id}</text>'
            )

    svg = (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}">\n'
        '<rect width="100%" height="100%" fill="white"/>\n'
        '<g stroke="#245b8f" stroke-width="1.8" fill="none">\n'
        + "\n".join(lines) +
        '\n</g>\n<g id="nodes">\n' + "\n".join(markers) +
        '\n</g>\n<g id="cell-labels">\n' + "\n".join(cell_labels) +
        '\n</g>\n</svg>\n'
    )
    path.write_text(svg, encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="OneFLOW structured .ofl file")
    parser.add_argument("--output", type=Path, default=Path("minimalRectangle.mesh"))
    parser.add_argument("--real-bytes", type=int, choices=(4, 8), default=8)
    parser.add_argument("--byte-order", choices=("little", "big"), default="little")
    args = parser.parse_args()
    endian = "<" if args.byte_order == "little" else ">"

    ni, nj, points = read_structured_grid(args.input, endian, args.real_bytes)
    cells, faces, boundary_faces = build_topology(ni, nj)
    output = args.output
    node_count = len(points)
    cell_count = len(cells)
    boundary_count = len(boundary_faces)

    write_json(output / "summary.json", {
        "schema": "oneflow.mesh-summary", "version": 1, "dimension": 2,
        "coordinateDimension": 3,
        "counts": {"nodes": node_count, "faces": len(faces), "cells": cell_count,
                   "boundaryFaces": boundary_count, "zones": 1},
        "boundingBox": coordinate_bounds(points),
        "source": str(args.input),
        "reconstructionLevel": "L1+L2",
    })
    write_json(output / "geometry.json", {
        "schema": "oneflow.mesh-geometry", "version": 1,
        "coordinateDimension": 3, "coordinateUnits": "unspecified",
        "idBase": 0, "points": points,
    })
    write_json(output / "topology.json", {
        "schema": "oneflow.mesh-topology", "version": 1,
        "idBase": 0, "cellTypes": {"quad": 9},
        "cells": cells, "faces": faces,
    })
    write_json(output / "boundary.json", {
        "schema": "oneflow.mesh-boundaries", "version": 1,
        "idBase": 0, "boundaryFaces": boundary_faces,
        "note": "The structured-grid source has no physical boundary names in this header.",
    })
    write_json(output / "manifest.json", {
        "schema": "oneflow.mesh-package", "version": 1,
        "sourceFormat": "oneflow-structured-binary",
        "idBase": 0,
        "files": ["summary.json", "geometry.json", "topology.json", "boundary.json",
                  "views/mesh.vtu", "views/mesh.svg"],
        "capabilities": {"geometry": True, "cellTopology": True, "faceTopology": True,
                         "physicalBoundaryNames": False, "partitionMaps": False},
    })
    write_vtu(output / "views" / "mesh.vtu", points, cells)
    write_svg(output / "views" / "mesh.svg", ni, nj, points)

    print(f"Exported mesh package: {output}")
    print(f"Dimension: 2; nodes: {node_count}; cells: {cell_count}; "
          f"faces: {len(faces)}; boundary faces: {boundary_count}")
    print(f"Views: {output / 'views' / 'mesh.svg'} and {output / 'views' / 'mesh.vtu'}")


if __name__ == "__main__":
    main()
