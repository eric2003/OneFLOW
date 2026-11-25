#!/usr/bin/env python3
"""Generate a deterministic, watertight 90-degree circular pipe-bend STL."""
from __future__ import annotations
import argparse
import hashlib
import json
import math
import struct
from collections import Counter
from pathlib import Path

REPOSITORY_ROOT = Path(__file__).resolve().parents[7]
DEFAULT_PARAMETERS = {
    "radius": 0.5,
    "bend_radius": 1.0,
    "inlet_length": 1.5,
    "outlet_length": 1.5,
    "angle_degrees": 90.0,
    "circumferential_divisions": 48,
    "bend_divisions": 24,
    "inlet_divisions": 12,
    "outlet_divisions": 12,
}


def _cross(a, b):
    return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])


def _sub(a, b):
    return tuple(a[i]-b[i] for i in range(3))


def _dot(a, b):
    return sum(a[i]*b[i] for i in range(3))


def _norm(a):
    return math.sqrt(_dot(a, a))


def build_surface(radius=0.5, bend_radius=1.0, inlet_length=1.5, outlet_length=1.5,
                  angle_degrees=90.0, circumferential_divisions=48,
                  bend_divisions=24, inlet_divisions=12, outlet_divisions=12):
    values = (radius, bend_radius, inlet_length, outlet_length, angle_degrees)
    if any(not math.isfinite(v) for v in values):
        raise ValueError("All geometric parameters must be finite.")
    if radius <= 0 or bend_radius <= radius or inlet_length <= 0 or outlet_length <= 0:
        raise ValueError("Require radius > 0, bend_radius > radius, and positive straight lengths.")
    if not 0 < angle_degrees <= 180:
        raise ValueError("angle_degrees must be in (0, 180].")
    if circumferential_divisions < 8 or bend_divisions < 2 or inlet_divisions < 1 or outlet_divisions < 1:
        raise ValueError("Use >=8 circumferential, >=2 bend, and >=1 straight divisions.")

    theta = math.radians(angle_degrees)
    rings = []
    for i in range(inlet_divisions + 1):
        s = -inlet_length + inlet_length * i / inlet_divisions
        rings.append(((s, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0)))
    for i in range(1, bend_divisions + 1):
        phi = theta * i / bend_divisions
        rings.append(((bend_radius * math.sin(phi), bend_radius * (1-math.cos(phi)), 0.0),
                      (math.cos(phi), math.sin(phi), 0.0),
                      (-math.sin(phi), math.cos(phi), 0.0)))
    end = (bend_radius * math.sin(theta), bend_radius * (1-math.cos(theta)), 0.0)
    tangent = (math.cos(theta), math.sin(theta), 0.0)
    normal = (-math.sin(theta), math.cos(theta), 0.0)
    for i in range(1, outlet_divisions + 1):
        s = outlet_length * i / outlet_divisions
        rings.append((tuple(end[k] + s*tangent[k] for k in range(3)), tangent, normal))

    n = circumferential_divisions
    vertices = []
    for center, _tangent, normal in rings:
        for j in range(n):
            alpha = 2 * math.pi * j / n
            vertices.append(tuple(center[k] + radius * (
                math.cos(alpha)*normal[k] + (math.sin(alpha) if k == 2 else 0.0)
            ) for k in range(3)))
    faces, labels = [], []
    for i in range(len(rings)-1):
        for j in range(n):
            a, b = i*n+j, i*n+(j+1)%n
            c, d = (i+1)*n+(j+1)%n, (i+1)*n+j
            faces.extend(((a, b, c), (a, c, d)))
            labels.extend(("wall", "wall"))
    inlet_center = len(vertices)
    vertices.append(rings[0][0])
    outlet_center = len(vertices)
    vertices.append(rings[-1][0])
    last = (len(rings)-1)*n
    for j in range(n):
        nxt = (j+1) % n
        faces.extend(((inlet_center, nxt, j), (outlet_center, last+j, last+nxt)))
        labels.extend(("inflow", "outflow"))

    signed = sum(_dot(vertices[a], _cross(vertices[b], vertices[c]))/6.0 for a,b,c in faces)
    if signed < 0:
        faces = [(a, c, b) for a,b,c in faces]
    return vertices, faces, labels


def inspect_surface(vertices, faces):
    uses, directed = Counter(), Counter()
    parent = list(range(len(vertices)))
    area = volume = 0.0

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for a,b,c in faces:
        for left, right in ((a, b), (b, c), (c, a)):
            rl, rr = find(left), find(right)
            if rl != rr:
                parent[rr] = rl
        pa,pb,pc = vertices[a],vertices[b],vertices[c]
        cross = _cross(_sub(pb,pa), _sub(pc,pa))
        area += 0.5*_norm(cross)
        volume += _dot(pa,_cross(pb,pc))/6.0
        for x,y in ((a,b),(b,c),(c,a)):
            uses[(min(x,y),max(x,y))] += 1
            directed[(x,y)] += 1
    bad = sum(count != 2 for count in uses.values())
    inconsistent = sum(1 for x,y in uses if directed[(x,y)] != 1 or directed[(y,x)] != 1)
    return {
        "vertex_count": len(vertices), "triangle_count": len(faces),
        "connected_component_count": len({find(i) for i in range(len(vertices))}),
        "boundary_or_nonmanifold_edge_count": bad,
        "inconsistently_oriented_edge_count": inconsistent,
        "surface_area": area, "signed_volume": volume,
        "watertight": bad == 0 and inconsistent == 0 and len({find(i) for i in range(len(vertices))}) == 1,
    }


def write_stl(path, vertices, faces):
    with path.open("wb") as stream:
        stream.write(b"OneFLOW deterministic 90-degree pipe bend".ljust(80, b"\0")[:80])
        stream.write(struct.pack("<I", len(faces)))
        for a,b,c in faces:
            pa,pb,pc = vertices[a],vertices[b],vertices[c]
            normal = _cross(_sub(pb,pa), _sub(pc,pa))
            length = _norm(normal)
            if length:
                normal = tuple(v/length for v in normal)
            stream.write(struct.pack("<12fH", *normal, *pa, *pb, *pc, 0))


def generate(output, *, units="m", **parameters):
    output = Path(output).expanduser().resolve()
    if output.suffix.lower() != ".stl":
        raise ValueError("Output path must use the .stl extension.")
    if output == REPOSITORY_ROOT or REPOSITORY_ROOT in output.parents:
        raise ValueError("Generated STL assets must be written outside the OneFLOW source tree.")
    if units not in {"m", "mm", "cm", "in"}:
        raise ValueError("Units must be one of: m, mm, cm, in.")
    parameters = {**DEFAULT_PARAMETERS, **parameters}
    vertices, faces, labels = build_surface(**parameters)
    metrics = inspect_surface(vertices, faces)
    if not metrics["watertight"] or metrics["signed_volume"] <= 0:
        raise ValueError(f"Generated surface failed topology/orientation checks: {metrics}")
    analytic_volume = math.pi*parameters["radius"]**2*(
        parameters["inlet_length"] + parameters["outlet_length"]
        + parameters["bend_radius"]*math.radians(parameters["angle_degrees"])
    )
    metrics["analytic_volume"] = analytic_volume
    metrics["relative_volume_error"] = abs(metrics["signed_volume"]-analytic_volume)/analytic_volume
    if metrics["relative_volume_error"] > 0.03:
        raise ValueError(f"STL volume differs excessively from analytic volume: {metrics}")
    output.parent.mkdir(parents=True, exist_ok=True)
    write_stl(output, vertices, faces)
    patch_indices = {}
    for index, label in enumerate(labels):
        patch_indices.setdefault(label, []).append(index)
    region_records = {}
    for name in ("wall", "inflow", "outflow"):
        region_path = output.parent / f"{output.stem}.{name}.stl"
        region_faces = [faces[index] for index in patch_indices[name]]
        write_stl(region_path, vertices, region_faces)
        region_records[name] = {
            "file": region_path.name,
            "sha256": hashlib.sha256(region_path.read_bytes()).hexdigest(),
            "triangle_count": len(region_faces),
        }
    region_manifest_path = output.parent / "region_manifest.json"
    region_manifest = {
        "schema": "oneflow.stl-surface-regions",
        "version": 1,
        "closed_surface_sha256": hashlib.sha256(output.read_bytes()).hexdigest(),
        "regions": region_records,
    }
    region_manifest_path.write_text(
        json.dumps(region_manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    definition = {
        "kind": "pipe_bend_90_degree", "units": units, "parameters": parameters,
    }
    definition_sha256 = hashlib.sha256(
        json.dumps(definition, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()
    generator_source = Path(__file__).resolve()
    generator_source_sha256 = hashlib.sha256(generator_source.read_bytes()).hexdigest()
    manifest = {
        "schema": "oneflow.stl-geometry-manifest", "version": 1,
        "asset_id": "geometry.pipe_bend.90_degree", "geometry": "pipe_bend_90_degree",
        "units": units, "parameters": parameters, "definition": definition,
        "definition_sha256": definition_sha256,
        "generator": {
            "name": "oneflow_pipe_bend_geometry", "version": "1.0",
            "source": generator_source.relative_to(REPOSITORY_ROOT).as_posix(),
            "source_sha256": generator_source_sha256,
        },
        "surface": metrics,
        "patch_triangle_indices": patch_indices,
        "region_manifest": region_manifest_path.name,
        "region_triangle_counts": {
            name: region_records[name]["triangle_count"] for name in ("wall", "inflow", "outflow")
        },
        "stl_file": output.name,
        "stl_sha256": hashlib.sha256(output.read_bytes()).hexdigest(),
        "intended_gmsh_physical_groups": {
            "surface": ["wall", "inflow", "outflow"], "volume": ["fluid_domain"]
        },
        "notes": "Inlet/outlet caps close the fluid-domain shell; they are boundary patches, not solid walls.",
    }
    sidecar = output.with_suffix(".geometry.json")
    sidecar.write_text(json.dumps(manifest, indent=2, sort_keys=True)+"\n", encoding="utf-8")
    return {
        "stl": output, "manifest": sidecar,
        "region_manifest_path": region_manifest_path,
        "region_files": {
            name: output.parent / region_records[name]["file"]
            for name in ("wall", "inflow", "outflow")
        },
        **manifest,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--radius", type=float, default=0.5)
    parser.add_argument("--bend-radius", type=float, default=1.0)
    parser.add_argument("--inlet-length", type=float, default=1.5)
    parser.add_argument("--outlet-length", type=float, default=1.5)
    parser.add_argument("--angle-degrees", type=float, default=90.0)
    parser.add_argument("--circumferential-divisions", type=int, default=48)
    parser.add_argument("--bend-divisions", type=int, default=24)
    parser.add_argument("--inlet-divisions", type=int, default=12)
    parser.add_argument("--outlet-divisions", type=int, default=12)
    parser.add_argument("--units", default="m")
    args = parser.parse_args()
    parameters = {key: value for key,value in vars(args).items() if key not in {"output", "units"}}
    try:
        report = generate(args.output, units=args.units, **parameters)
    except (OSError, ValueError) as error:
        parser.exit(2, f"Pipe-bend generation failed: {error}\n")
    print(json.dumps(report, indent=2, default=str))


if __name__ == "__main__":
    main()
