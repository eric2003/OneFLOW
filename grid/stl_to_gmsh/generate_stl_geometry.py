#!/usr/bin/env python3
"""Deterministic parameterized STL surface generation for Mesh Generation."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import struct
from collections import Counter
from pathlib import Path

GENERATOR_NAME = "oneflow_stl_geometry"
GENERATOR_VERSION = "1.0"


def _orient_outward(vertices: list[tuple[float, float, float]], faces: list[tuple[int, int, int]]) -> list[tuple[int, int, int]]:
    oriented = []
    for a, b, c in faces:
        pa, pb, pc = vertices[a], vertices[b], vertices[c]
        ab = tuple(pb[i] - pa[i] for i in range(3))
        ac = tuple(pc[i] - pa[i] for i in range(3))
        normal = (
            ab[1] * ac[2] - ab[2] * ac[1],
            ab[2] * ac[0] - ab[0] * ac[2],
            ab[0] * ac[1] - ab[1] * ac[0],
        )
        center = tuple((pa[i] + pb[i] + pc[i]) / 3.0 for i in range(3))
        if sum(normal[i] * center[i] for i in range(3)) < 0:
            oriented.append((a, c, b))
        else:
            oriented.append((a, b, c))
    return oriented


def _icosphere(subdivisions: int, radius: float) -> tuple[list[tuple[float, float, float]], list[tuple[int, int, int]]]:
    if not 0 <= subdivisions <= 6:
        raise ValueError("Sphere subdivision level must be between 0 and 6.")
    if not math.isfinite(radius) or radius <= 0:
        raise ValueError("Sphere radius must be positive and finite.")
    phi = (1.0 + math.sqrt(5.0)) / 2.0
    raw = [
        (-1, phi, 0), (1, phi, 0), (-1, -phi, 0), (1, -phi, 0),
        (0, -1, phi), (0, 1, phi), (0, -1, -phi), (0, 1, -phi),
        (phi, 0, -1), (phi, 0, 1), (-phi, 0, -1), (-phi, 0, 1),
    ]
    vertices = []
    for point in raw:
        length = math.sqrt(sum(value * value for value in point))
        vertices.append(tuple(radius * value / length for value in point))
    faces = [
        (0, 11, 5), (0, 5, 1), (0, 1, 7), (0, 7, 10), (0, 10, 11),
        (1, 5, 9), (5, 11, 4), (11, 10, 2), (10, 7, 6), (7, 1, 8),
        (3, 9, 4), (3, 4, 2), (3, 2, 6), (3, 6, 8), (3, 8, 9),
        (4, 9, 5), (2, 4, 11), (6, 2, 10), (8, 6, 7), (9, 8, 1),
    ]
    faces = _orient_outward(vertices, faces)
    for _ in range(subdivisions):
        midpoint_cache: dict[tuple[int, int], int] = {}

        def midpoint(left: int, right: int) -> int:
            key = (min(left, right), max(left, right))
            if key not in midpoint_cache:
                a, b = vertices[left], vertices[right]
                point = tuple((a[i] + b[i]) * 0.5 for i in range(3))
                length = math.sqrt(sum(value * value for value in point))
                midpoint_cache[key] = len(vertices)
                vertices.append(tuple(radius * value / length for value in point))
            return midpoint_cache[key]

        refined = []
        for a, b, c in faces:
            ab, bc, ca = midpoint(a, b), midpoint(b, c), midpoint(c, a)
            refined.extend(((a, ab, ca), (b, bc, ab), (c, ca, bc), (ab, bc, ca)))
        faces = _orient_outward(vertices, refined)
    return vertices, faces


def _cylinder(radius: float, length: float, circumferential_divisions: int, axial_divisions: int):
    if not math.isfinite(radius) or radius <= 0 or not math.isfinite(length) or length <= 0:
        raise ValueError("Cylinder radius and length must be positive and finite.")
    if circumferential_divisions < 8 or axial_divisions < 1:
        raise ValueError("Cylinder needs at least 8 circumferential divisions and 1 axial interval.")
    vertices: list[tuple[float, float, float]] = []
    for i in range(axial_divisions + 1):
        x = -0.5 * length + length * i / axial_divisions
        for j in range(circumferential_divisions):
            angle = 2.0 * math.pi * j / circumferential_divisions
            vertices.append((x, radius * math.cos(angle), radius * math.sin(angle)))
    faces: list[tuple[int, int, int]] = []
    labels: list[str] = []
    n = circumferential_divisions
    for i in range(axial_divisions):
        for j in range(n):
            a = i * n + j
            b = i * n + (j + 1) % n
            c = (i + 1) * n + (j + 1) % n
            d = (i + 1) * n + j
            faces.extend(((a, b, c), (a, c, d)))
            labels.extend(("wall", "wall"))
    inlet_center = len(vertices)
    vertices.append((-0.5 * length, 0.0, 0.0))
    outlet_center = len(vertices)
    vertices.append((0.5 * length, 0.0, 0.0))
    for j in range(n):
        nxt = (j + 1) % n
        faces.append((inlet_center, nxt, j))
        labels.append("inflow")
        outlet_base = axial_divisions * n
        faces.append((outlet_center, outlet_base + j, outlet_base + nxt))
        labels.append("outflow")
    # The construction order is deliberate; verify winding rather than silently
    # reorienting patches because the cap normals carry boundary semantics.
    return vertices, faces, labels


def _surface_metrics(vertices: list[tuple[float, float, float]], faces: list[tuple[int, int, int]]) -> dict:
    edge_uses: Counter[tuple[int, int]] = Counter()
    directed: Counter[tuple[int, int]] = Counter()
    area = 0.0
    signed_volume = 0.0
    parent = list(range(len(vertices)))

    def find(item: int) -> int:
        while parent[item] != item:
            parent[item] = parent[parent[item]]
            item = parent[item]
        return item

    def union(left: int, right: int) -> None:
        left_root, right_root = find(left), find(right)
        if left_root != right_root:
            parent[right_root] = left_root

    for a, b, c in faces:
        union(a, b)
        union(b, c)
        pa, pb, pc = vertices[a], vertices[b], vertices[c]
        ab = tuple(pb[i] - pa[i] for i in range(3))
        ac = tuple(pc[i] - pa[i] for i in range(3))
        cross = (
            ab[1] * ac[2] - ab[2] * ac[1],
            ab[2] * ac[0] - ab[0] * ac[2],
            ab[0] * ac[1] - ab[1] * ac[0],
        )
        area += 0.5 * math.sqrt(sum(value * value for value in cross))
        signed_volume += (
            pa[0] * (pb[1] * pc[2] - pb[2] * pc[1])
            + pa[1] * (pb[2] * pc[0] - pb[0] * pc[2])
            + pa[2] * (pb[0] * pc[1] - pb[1] * pc[0])
        ) / 6.0
        for left, right in ((a, b), (b, c), (c, a)):
            key = (min(left, right), max(left, right))
            edge_uses[key] += 1
            directed[(left, right)] += 1
    bad_edges = sum(count != 2 for count in edge_uses.values())
    component_count = len({find(index) for index in range(len(vertices))})
    inconsistent = sum(
        1 for left, right in edge_uses
        if directed[(left, right)] + directed[(right, left)] != 2
        or directed[(left, right)] != 1
        or directed[(right, left)] != 1
    )
    bounds = {
        axis: [min(point[i] for point in vertices), max(point[i] for point in vertices)]
        for i, axis in enumerate(("x", "y", "z"))
    }
    return {
        "vertex_count": len(vertices),
        "triangle_count": len(faces),
        "connected_component_count": component_count,
        "surface_area": area,
        "signed_volume": signed_volume,
        "bounds": bounds,
        "boundary_or_nonmanifold_edge_count": bad_edges,
        "inconsistently_oriented_edge_count": inconsistent,
        "watertight": bad_edges == 0 and inconsistent == 0,
    }


def _write_binary_stl(path: Path, vertices: list[tuple[float, float, float]], faces: list[tuple[int, int, int]]) -> None:
    with path.open("wb") as stream:
        header = b"OneFLOW deterministic STL geometry generator"
        stream.write(header.ljust(80, b"\0")[:80])
        stream.write(struct.pack("<I", len(faces)))
        for a, b, c in faces:
            pa, pb, pc = vertices[a], vertices[b], vertices[c]
            ab = tuple(pb[i] - pa[i] for i in range(3))
            ac = tuple(pc[i] - pa[i] for i in range(3))
            normal = (
                ab[1] * ac[2] - ab[2] * ac[1],
                ab[2] * ac[0] - ab[0] * ac[2],
                ab[0] * ac[1] - ab[1] * ac[0],
            )
            norm = math.sqrt(sum(value * value for value in normal))
            if norm:
                normal = tuple(value / norm for value in normal)
            stream.write(struct.pack("<12fH", *normal, *vertices[a], *vertices[b], *vertices[c], 0))


def generate_geometry(
    geometry: str,
    output_stl: Path,
    *,
    radius: float = 1.0,
    subdivisions: int = 3,
    semi_axes: tuple[float, float, float] = (3.0, 0.5, 0.5),
    length: float = 4.0,
    circumferential_divisions: int = 48,
    axial_divisions: int = 16,
    units: str = "m",
) -> dict:
    """Generate a closed STL and checksum-bound metadata manifest."""
    output_stl = Path(output_stl).expanduser().resolve()
    if output_stl.suffix.lower() != ".stl":
        raise ValueError("STL output path must use the .stl extension.")
    if not units.strip():
        raise ValueError("A geometry unit must be declared.")
    if geometry in {"sphere", "prolate_spheroid"}:
        vertices, faces = _icosphere(subdivisions, radius if geometry == "sphere" else 1.0)
        if geometry == "prolate_spheroid":
            if len(semi_axes) != 3 or any(not math.isfinite(v) or v <= 0 for v in semi_axes):
                raise ValueError("All spheroid semi-axes must be positive and finite.")
            vertices = [
                (point[0] * semi_axes[0], point[1] * semi_axes[1], point[2] * semi_axes[2])
                for point in vertices
            ]
        patch_labels = ["wall"] * len(faces)
        parameters = {
            **({"radius": radius} if geometry == "sphere" else {"semi_axes": list(semi_axes)}),
            "subdivisions": subdivisions,
        }
    elif geometry == "cylinder_internal_flow":
        vertices, faces, patch_labels = _cylinder(
            radius, length, circumferential_divisions, axial_divisions
        )
        parameters = {
            "radius": radius, "length": length,
            "circumferential_divisions": circumferential_divisions,
            "axial_divisions": axial_divisions,
        }
    else:
        raise ValueError(f"Unsupported geometry type: {geometry}")

    metrics = _surface_metrics(vertices, faces)
    if not metrics["watertight"] or metrics["signed_volume"] <= 0:
        raise ValueError(
            "Generated surface failed closure/orientation checks: "
            f"watertight={metrics['watertight']}, signed_volume={metrics['signed_volume']}"
        )
    output_stl.parent.mkdir(parents=True, exist_ok=True)
    _write_binary_stl(output_stl, vertices, faces)
    digest = hashlib.sha256(output_stl.read_bytes()).hexdigest()
    patch_ranges: dict[str, list[int]] = {}
    for index, label in enumerate(patch_labels):
        patch_ranges.setdefault(label, []).append(index)
    manifest = {
        "schema": "oneflow.stl-geometry-manifest",
        "version": 1,
        "generator": GENERATOR_NAME,
        "generator_version": GENERATOR_VERSION,
        "geometry": geometry,
        "units": units,
        "parameters": parameters,
        "surface": metrics,
        "patch_triangle_indices": patch_ranges,
        "stl_file": output_stl.name,
        "stl_sha256": digest,
        "physical_groups": (
            {"surface": ["wall", "inflow", "outflow"], "volume": ["fluid_domain"]}
            if geometry == "cylinder_internal_flow"
            else {"surface": ["wall"], "volume": ["fluid_domain"]}
        ),
        "notes": (
            "STL facets are geometric closure patches. Inflow and outflow are not solid walls."
            if geometry == "cylinder_internal_flow"
            else "The STL is a piecewise-planar approximation of the analytic surface."
        ),
    }
    manifest_path = output_stl.with_suffix(".geometry.json")
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    return {"stl": output_stl, "manifest": manifest_path, **manifest}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--geometry", required=True, choices=("sphere", "prolate_spheroid", "cylinder_internal_flow"))
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--radius", type=float, default=1.0)
    parser.add_argument("--subdivisions", type=int, default=3)
    parser.add_argument("--semi-axis-a", type=float, default=3.0)
    parser.add_argument("--semi-axis-b", type=float, default=0.5)
    parser.add_argument("--semi-axis-c", type=float, default=0.5)
    parser.add_argument("--length", type=float, default=4.0)
    parser.add_argument("--circumferential-divisions", type=int, default=48)
    parser.add_argument("--axial-divisions", type=int, default=16)
    parser.add_argument("--units", default="m")
    args = parser.parse_args()
    try:
        report = generate_geometry(
            args.geometry, args.output, radius=args.radius,
            subdivisions=args.subdivisions,
            semi_axes=(args.semi_axis_a, args.semi_axis_b, args.semi_axis_c),
            length=args.length,
            circumferential_divisions=args.circumferential_divisions,
            axial_divisions=args.axial_divisions, units=args.units,
        )
    except (OSError, ValueError) as error:
        parser.exit(2, f"STL geometry generation failed: {error}\n")
    print(json.dumps(report, indent=2, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
