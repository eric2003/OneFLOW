#!/usr/bin/env python3
"""Generate the OECD/NEA-Vattenfall T-junction surface and named STL regions.

The Boolean fluid volume is constructed by Gmsh/OpenCASCADE. The resulting
MSH2 surface triangles are exported as a closed STL plus checksum-bound
wall/inflow/outflow compatibility patches plus separate main/branch inlet subregions for the shared OneFLOW region-aware STL mesher.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import shutil
import struct
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[7]
CASE_DIRECTORY = Path(__file__).resolve().parents[2] / "automatic_tetrahedra"
_BENCHMARK = json.loads((CASE_DIRECTORY / "benchmark.json").read_text(encoding="utf-8"))
_GEOMETRY = _BENCHMARK["geometry"]
DEFAULTS = {
    "main_diameter": _GEOMETRY["main_pipe"]["inner_diameter"],
    "branch_diameter": _GEOMETRY["branch_pipe"]["inner_diameter"],
    "main_inlet_length": _GEOMETRY["main_pipe"]["inlet_length"],
    "main_outlet_length": _GEOMETRY["main_pipe"]["outlet_length"],
    "branch_inlet_length": _GEOMETRY["branch_pipe"]["inlet_length"],
    "mesh_size": 0.030,
}
GROUPS = {"wall": 1, "main_inflow": 2, "branch_inflow": 3, "outflow": 4}


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _geo_text(p: dict) -> str:
    rm = p["main_diameter"] / 2.0
    rb = p["branch_diameter"] / 2.0
    li, lo, lb = p["main_inlet_length"], p["main_outlet_length"], p["branch_inlet_length"]
    # OCC Boolean operations introduce small bounding-box tolerances. Keep
    # the selector slab wider than the kernel tolerance and scale it mildly
    # with the requested mesh size.
    eps = max(1.0e-5, p["mesh_size"] * 1.0e-3)
    return f"""SetFactory("OpenCASCADE");
Geometry.OCCBooleanPreserveNumbering = 0;
Mesh.MshFileVersion = 2.2;
Mesh.MeshSizeMin = {p["mesh_size"]:.12g};
Mesh.MeshSizeMax = {p["mesh_size"]:.12g};
Mesh.MeshSizeFromCurvature = 16;
Mesh.MeshSizeExtendFromBoundary = 1;
Cylinder(1) = {{{-li:.12g}, 0, 0, {li + lo:.12g}, 0, 0, {rm:.12g}}};
Cylinder(2) = {{0, 0, 0, 0, 0, {lb:.12g}, {rb:.12g}}};
fluid[] = BooleanUnion{{ Volume{{1}}; Delete; }}{{ Volume{{2}}; Delete; }};
surfaces[] = Boundary{{ Volume{{fluid[]}}; }};
inlet[] = Surface In BoundingBox{{ {-li-eps:.12g}, {-rm-eps:.12g}, {-rm-eps:.12g}, {-li+eps:.12g}, {rm+eps:.12g}, {rm+eps:.12g} }};
outlet[] = Surface In BoundingBox{{ {lo-eps:.12g}, {-rm-eps:.12g}, {-rm-eps:.12g}, {lo+eps:.12g}, {rm+eps:.12g}, {rm+eps:.12g} }};
branch_inlet[] = Surface In BoundingBox{{ {-rb-eps:.12g}, {-rb-eps:.12g}, {lb-eps:.12g}, {rb+eps:.12g}, {rb+eps:.12g}, {lb+eps:.12g} }};
wall[] = surfaces[];
wall[] -= inlet[];
wall[] -= outlet[];
wall[] -= branch_inlet[];
Physical Surface("wall", 1) = {{wall[]}};
Physical Surface("main_inflow", 2) = {{inlet[]}};
Physical Surface("branch_inflow", 3) = {{branch_inlet[]}};
Physical Surface("outflow", 4) = {{outlet[]}};
Physical Volume("fluid_domain", 1) = {{fluid[]}};
"""


def _validate_pipe_tee_surface_geometry(nodes: dict, groups: dict, parameters: dict) -> dict:
    """Check the three orthogonal circular ports and cylindrical wall of the T-junction."""
    required = ("wall", "main_inflow", "branch_inflow", "outflow")
    if any(not groups.get(name) for name in required):
        raise ValueError("T-junction geometry is missing one or more physical surface groups.")
    rm = parameters["main_diameter"] / 2.0
    rb = parameters["branch_diameter"] / 2.0
    li, lo, lb = parameters["main_inlet_length"], parameters["main_outlet_length"], parameters["branch_inlet_length"]
    mesh_size = parameters["mesh_size"]
    plane_tol = max(1.0e-8, mesh_size * 1.0e-4)
    radial_tol = max(1.0e-5, mesh_size * 1.0e-3)
    specs = (
        ("main_inflow", 0, -li, 1, 2, rm, (-li, 0.0, 0.0), "+X"),
        ("outflow", 0, lo, 1, 2, rm, (lo, 0.0, 0.0), "+X"),
        ("branch_inflow", 2, lb, 0, 1, rb, (0.0, 0.0, lb), "+Z"),
    )
    ports = {}
    for name, axis_index, plane, radial_a, radial_b, radius, center, axis in specs:
        area = 0.0
        centroid_sum = [0.0, 0.0, 0.0]
        max_radius = 0.0
        for face in groups[name]:
            pts = [nodes[i] for i in face]
            if any(abs(p[axis_index] - plane) > plane_tol for p in pts):
                raise ValueError(f"T-junction port '{name}' is not on its expected {axis} end plane.")
            a, b, c = pts
            ab = [b[i] - a[i] for i in range(3)]
            ac = [c[i] - a[i] for i in range(3)]
            cross = (ab[1]*ac[2]-ab[2]*ac[1], ab[2]*ac[0]-ab[0]*ac[2], ab[0]*ac[1]-ab[1]*ac[0])
            tri_area = 0.5 * math.sqrt(sum(v*v for v in cross))
            area += tri_area
            for i in range(3):
                centroid_sum[i] += tri_area * sum(p[i] for p in pts) / 3.0
            for p in pts:
                r = math.hypot(p[radial_a]-center[radial_a], p[radial_b]-center[radial_b])
                max_radius = max(max_radius, r)
        if area <= 0.0:
            raise ValueError(f"T-junction port '{name}' has zero surface area.")
        centroid = [v / area for v in centroid_sum]
        center_offset = math.hypot(centroid[radial_a]-center[radial_a], centroid[radial_b]-center[radial_b])
        expected_area = math.pi * radius * radius
        area_error = abs(area-expected_area) / expected_area
        if area_error > 0.20:
            raise ValueError(f"T-junction port '{name}' area differs from its circular design area by {area_error:.1%}.")
        if center_offset > radius * 0.15:
            raise ValueError(f"T-junction port '{name}' is not centered on its design axis.")
        if not radius * 0.999 <= max_radius <= radius * 1.001:
            raise ValueError(f"T-junction port '{name}' boundary radius does not match its specified radius.")
        ports[name] = {
            "axis": axis, "plane_coordinate_m": plane, "measured_area_m2": area,
            "expected_area_m2": expected_area, "relative_area_error": area_error,
            "area_centroid_m": centroid, "maximum_radius_m": max_radius,
            "triangle_count": len(groups[name]),
        }
    wall_ids = {i for face in groups["wall"] for i in face}
    if not wall_ids:
        raise ValueError("T-junction wall Physical Surface contains no mesh nodes.")
    wall_surface_counts = {"main_pipe": 0, "branch_pipe": 0}
    checked_wall_faces = 0
    for face in groups["wall"]:
        pts = [nodes[i] for i in face]
        main_residuals = [abs(math.hypot(p[1], p[2]) - rm) for p in pts]
        branch_residuals = [abs(math.hypot(p[0], p[1]) - rb) for p in pts]
        on_main = max(main_residuals) <= radial_tol
        on_branch = max(branch_residuals) <= radial_tol
        if on_main == on_branch:
            raise ValueError(
                "A wall triangle must belong to exactly one analytic cylinder surface; "
                "it is either off-surface or ambiguously on the cylinder intersection."
            )
        centroid = tuple(sum(p[axis] for p in pts) / 3.0 for axis in range(3))
        x, y, z = centroid
        if on_main:
            inside_other_bore = (
                0.0 + radial_tol < z < lb - radial_tol
                and x * x + y * y < (rb - radial_tol) ** 2
            )
            if inside_other_bore:
                raise ValueError(
                    "A main-pipe wall triangle lies inside the branch bore; "
                    "BooleanUnion must remove this internal interface from wall."
                )
            wall_surface_counts["main_pipe"] += 1
        else:
            inside_other_bore = (
                -li + radial_tol < x < lo - radial_tol
                and y * y + z * z < (rm - radial_tol) ** 2
            )
            if inside_other_bore:
                raise ValueError(
                    "A branch-pipe wall triangle lies inside the main bore; "
                    "BooleanUnion must remove this internal interface from wall."
                )
            wall_surface_counts["branch_pipe"] += 1
        checked_wall_faces += 1
    if not wall_surface_counts["main_pipe"] or not wall_surface_counts["branch_pipe"]:
        raise ValueError("T-junction wall must contain exposed faces from both pipe cylinders.")
    return {
        "status": "PASS", "geometry": "orthogonal_three_port_t_junction", "units": "m",
        "main_axis": "+X", "branch_axis": "+Z", "axis_angle_degrees": 90.0,
        "port_count": len(ports), "ports": ports, "wall_node_count": len(wall_ids),
        "wall_triangle_count": checked_wall_faces,
        "wall_surface_triangle_counts": wall_surface_counts,
        "wall_radial_tolerance_m": radial_tol,
        "internal_overlap_wall_triangles": 0,
    }


def _parse_msh2(path: Path) -> tuple[dict[int, tuple[float, float, float]], dict[str, list[tuple[int, int, int]]]]:
    lines = path.read_text(encoding="utf-8", errors="strict").splitlines()
    physical_names: dict[tuple[int, int], str] = {}
    nodes: dict[int, tuple[float, float, float]] = {}
    triangles: dict[str, list[tuple[int, int, int]]] = {name: [] for name in GROUPS}
    i = 0
    while i < len(lines):
        marker = lines[i].strip()
        if marker == "$PhysicalNames":
            count = int(lines[i + 1].strip())
            for line in lines[i + 2:i + 2 + count]:
                dim_text, tag_text, name_text = line.split(maxsplit=2)
                physical_names[(int(dim_text), int(tag_text))] = name_text.strip().strip('"')
            i += count + 3
            continue
        if marker == "$Nodes":
            count = int(lines[i + 1].strip())
            for line in lines[i + 2:i + 2 + count]:
                fields = line.split()
                nodes[int(fields[0])] = tuple(map(float, fields[1:4]))
            i += count + 3
            continue
        if marker == "$Elements":
            count = int(lines[i + 1].strip())
            for line in lines[i + 2:i + 2 + count]:
                fields = list(map(int, line.split()))
                element_type, tag_count = fields[1], fields[2]
                if element_type != 2 or tag_count < 1:
                    continue
                name = physical_names.get((2, fields[3]))
                if name in triangles:
                    node_ids = fields[3 + tag_count:]
                    if len(node_ids) != 3:
                        raise ValueError("Unexpected non-triangular surface element in MSH2 output.")
                    triangles[name].append(tuple(node_ids))
            i += count + 3
            continue
        i += 1
    if not nodes:
        raise ValueError(f"No nodes were found in Gmsh surface mesh: {path}")
    if any(not triangles[name] for name in GROUPS):
        counts = {name: len(items) for name, items in triangles.items()}
        raise ValueError(f"One or more T-junction Physical Surface groups are empty: {counts}")
    return nodes, triangles


def _write_stl(path: Path, nodes: dict, faces: list[tuple[int, int, int]], label: str) -> None:
    with path.open("wb") as stream:
        stream.write((f"OneFLOW OECD Vattenfall T-junction {label}").encode("ascii").ljust(80, b"\0")[:80])
        stream.write(struct.pack("<I", len(faces)))
        for face in faces:
            a, b, c = (nodes[index] for index in face)
            ab = (b[0]-a[0], b[1]-a[1], b[2]-a[2])
            ac = (c[0]-a[0], c[1]-a[1], c[2]-a[2])
            normal = (ab[1]*ac[2]-ab[2]*ac[1], ab[2]*ac[0]-ab[0]*ac[2], ab[0]*ac[1]-ab[1]*ac[0])
            length = math.sqrt(sum(value*value for value in normal))
            if length:
                normal = tuple(value/length for value in normal)
            stream.write(struct.pack("<12fH", *normal, *a, *b, *c, 0))


def export_surface(msh: Path, output: Path, *, generator_kind: str = "python",
                   generator_source: Path | None = None, mesh_size: float | None = None,
                   geometry_parameters: dict | None = None,
                   operating_mode: str = "benchmark") -> dict:
    output = output.expanduser().resolve()
    if output.suffix.lower() != ".stl":
        raise ValueError("T-junction output must use the .stl extension.")
    if output == ROOT or ROOT in output.parents:
        raise ValueError("Generated STL assets must stay outside the OneFLOW source tree.")
    output.parent.mkdir(parents=True, exist_ok=True)
    msh = Path(msh).expanduser().resolve()
    nodes, groups = _parse_msh2(msh)
    geometry_validation = (_validate_pipe_tee_surface_geometry(nodes, groups, geometry_parameters)
                           if geometry_parameters is not None else None)
    if operating_mode not in {"benchmark", "main_inlets", "branch_inlet"}:
        raise ValueError(f"Unsupported T-junction operating mode: {operating_mode}")
    port_faces = {
        "main_left": groups["main_inflow"],
        "main_right": groups["outflow"],
        "branch": groups["branch_inflow"],
    }
    if operating_mode == "main_inlets":
        final_regions = {
            "wall": groups["wall"],
            "inflow": port_faces["main_left"] + port_faces["main_right"],
            "outflow": port_faces["branch"],
        }
        named_subregions = {
            "main_left_inflow": ("inflow", port_faces["main_left"]),
            "main_right_inflow": ("inflow", port_faces["main_right"]),
            "branch_outflow": ("outflow", port_faces["branch"]),
        }
    elif operating_mode == "branch_inlet":
        final_regions = {
            "wall": groups["wall"],
            "inflow": port_faces["branch"],
            "outflow": port_faces["main_left"] + port_faces["main_right"],
        }
        named_subregions = {
            "branch_inflow": ("inflow", port_faces["branch"]),
            "main_left_outflow": ("outflow", port_faces["main_left"]),
            "main_right_outflow": ("outflow", port_faces["main_right"]),
        }
    else:
        final_regions = {
            "wall": groups["wall"],
            "inflow": port_faces["main_left"] + port_faces["branch"],
            "outflow": port_faces["main_right"],
        }
        named_subregions = {
            "main_inflow": ("inflow", port_faces["main_left"]),
            "branch_inflow": ("inflow", port_faces["branch"]),
        }
    all_faces = final_regions["wall"] + final_regions["inflow"] + final_regions["outflow"]
    _write_stl(output, nodes, all_faces, "closed-fluid-domain")
    records = {}
    for name, faces in final_regions.items():
        patch = output.parent / f"{output.stem}.{name}.stl"
        _write_stl(patch, nodes, faces, name)
        records[name] = {"file": patch.name, "sha256": _sha(patch), "triangle_count": len(faces)}
    subregions = {}
    for name, (parent_region, faces) in named_subregions.items():
        patch = output.parent / f"{output.stem}.{name}.stl"
        _write_stl(patch, nodes, faces, name)
        subregions[name] = {
            "file": patch.name, "sha256": _sha(patch),
            "triangle_count": len(faces), "parent_region": parent_region,
        }
    region_path = output.parent / "region_manifest.json"
    region_manifest = {
        "schema": "oneflow.stl-surface-regions",
        "version": 1,
        "closed_surface_sha256": _sha(output),
        "regions": records,
        "subregions": subregions,
    }
    region_path.write_text(json.dumps(region_manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    parameters = dict(DEFAULTS)
    parameters.pop("mesh_size")
    definition = {
        "kind": "oecd_nea_vattenfall_t_junction",
        "units": "m",
        "parameters": parameters,
        "coordinate_system": "main pipe axis +X; branch inlet axis +Z; junction origin at pipe centreline intersection",
    }
    definition_hash = hashlib.sha256(json.dumps(definition, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    source_path = (generator_source or Path(__file__)).expanduser().resolve()
    manifest = {
        "schema": "oneflow.stl-geometry-manifest",
        "version": 1,
        "asset_id": "geometry.pipe_tee.oecd_nea_vattenfall",
        "geometry": "pipe_tee_oecd_nea_vattenfall",
        "units": "m",
        "parameters": parameters,
        "operating_mode": operating_mode,
        "definition": definition,
        "definition_sha256": definition_hash,
        "generator": {
            "name": "oneflow_pipe_tee_geometry",
            "version": "1.0",
            "kind": generator_kind,
            "source": source_path.relative_to(ROOT).as_posix(),
            "source_sha256": _sha(source_path),
        },
        "mesh_generation": {
            "surface_mesh_size_m": mesh_size if mesh_size is not None else DEFAULTS["mesh_size"],
            "surface_mesh_file": msh.name,
            "surface_mesh_sha256": _sha(msh),
        },
        "geometry_validation": geometry_validation,
        "surface": {
            "node_count": len(nodes),
            "triangle_count": len(all_faces),
            "region_triangle_counts": {name: len(faces) for name, faces in final_regions.items()},
            "port_subregion_triangle_counts": {name: len(faces) for name, (_parent, faces) in named_subregions.items()},
        },
        "region_manifest": region_path.name,
        "stl_file": output.name,
        "stl_sha256": _sha(output),
        "intended_gmsh_physical_groups": {
            "surface": ["wall", "inflow", "outflow"],
            "surface_subregions": sorted(named_subregions),
            "direct_occ_volume_mesh_surface_groups": ["wall", "main_inflow", "branch_inflow", "outflow"],
            "volume": ["fluid_domain"],
        },
        "notes": "Fluid domain is the Boolean union of the main and branch cylindrical bores. The STL region manifest maps the selected operating mode to wall/inflow/outflow and retains each physical port as a checksum-bound subregion. The generic region-aware STL route reconstructs the fluid volume from the exported STL; the benchmark runner also retains a separate direct OpenCASCADE volume-mesh route. The benchmark junction edge is reported as very slightly rounded (radius < 0.2 mm, not measured); the baseline geometry uses a sharp intersection and does not invent a fillet radius.",
    }
    sidecar = output.with_suffix(".geometry.json")
    sidecar.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {"stl": output, "manifest": sidecar, "region_manifest": region_path, "surface": manifest["surface"]}


def generate(output: Path, gmsh: str, mesh_size: float = DEFAULTS["mesh_size"],
             operating_mode: str = "benchmark") -> dict:
    if not math.isfinite(mesh_size) or mesh_size <= 0:
        raise ValueError("mesh_size must be positive and finite.")
    executable = shutil.which(gmsh) or (str(Path(gmsh).resolve()) if Path(gmsh).is_file() else None)
    if not executable:
        raise FileNotFoundError(f"Gmsh executable was not found: {gmsh}")
    output = output.expanduser().resolve()
    if output == ROOT or ROOT in output.parents:
        raise ValueError("Generated STL assets must stay outside the OneFLOW source tree.")
    output.parent.mkdir(parents=True, exist_ok=True)
    parameters = {**DEFAULTS, "mesh_size": mesh_size}
    geo = output.with_suffix(".geo")
    msh = output.with_suffix(".surface.msh")
    geo.write_text(_geo_text(parameters), encoding="ascii")
    subprocess.run(
        [executable, str(geo), "-2", "-format", "msh2", "-o", str(msh), "-v", "2"],
        check=True, cwd=output.parent, capture_output=True, text=True,
    )
    return export_surface(msh, output, mesh_size=mesh_size, geometry_parameters=parameters,
                          operating_mode=operating_mode)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--gmsh")
    parser.add_argument("--mesh-size", type=float, default=DEFAULTS["mesh_size"])
    parser.add_argument("--operating-mode", choices=("benchmark", "main_inlets", "branch_inlet"), default="benchmark")
    parser.add_argument("--extract-msh", type=Path)
    parser.add_argument("--generator-kind", choices=("python", "cpp"), default="python")
    parser.add_argument("--generator-source", type=Path)
    parser.add_argument("--validate-pipe-tee", action="store_true")
    parser.add_argument("--main-diameter", type=float, default=DEFAULTS["main_diameter"])
    parser.add_argument("--branch-diameter", type=float, default=DEFAULTS["branch_diameter"])
    parser.add_argument("--main-inlet-length", type=float, default=DEFAULTS["main_inlet_length"])
    parser.add_argument("--main-outlet-length", type=float, default=DEFAULTS["main_outlet_length"])
    parser.add_argument("--branch-inlet-length", type=float, default=DEFAULTS["branch_inlet_length"])
    args = parser.parse_args()
    try:
        geometry_parameters = ({"main_diameter": args.main_diameter, "branch_diameter": args.branch_diameter,
                                "main_inlet_length": args.main_inlet_length, "main_outlet_length": args.main_outlet_length,
                                "branch_inlet_length": args.branch_inlet_length, "mesh_size": args.mesh_size}
                               if args.validate_pipe_tee else None)
        result = export_surface(args.extract_msh, args.output, generator_kind=args.generator_kind,
                                generator_source=args.generator_source, mesh_size=args.mesh_size,
                                geometry_parameters=geometry_parameters,
                                operating_mode=args.operating_mode) if args.extract_msh else generate(
                                    args.output, args.gmsh, args.mesh_size, args.operating_mode
                                )
    except (OSError, ValueError, subprocess.CalledProcessError) as error:
        parser.exit(2, f"T-junction geometry generation failed: {error}\n")
    print(json.dumps({key: str(value) if isinstance(value, Path) else value for key, value in result.items()}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
