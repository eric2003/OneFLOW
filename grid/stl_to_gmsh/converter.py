"""Convert STL triangle surfaces into Gmsh geometry and inspect MSH output."""
from __future__ import annotations
import argparse, json, math, re, subprocess
from collections import Counter, defaultdict
from pathlib import Path

from grid.stl.reader import Point, Triangle, StlReadError, read_stl as _read_stl
from grid.stl.surface import StlGeometryError, _weld_vertices, extract_boundary

Stl2DError = StlGeometryError

def read_stl(path: Path) -> list[Triangle]:
    """Preserve the converter's historical exception type for STL read errors."""
    try:
        return _read_stl(path)
    except StlReadError as exc:
        raise Stl2DError(str(exc)) from exc

def _read_sidecar(path: Path, loop: list[Point], tol: float) -> tuple[dict, dict[tuple[int,int],str]]:
    config = json.loads(path.read_text(encoding="utf-8"))
    if "domain" in config:
        domain = config["domain"]
        if not isinstance(domain, dict):
            raise Stl2DError("domain must be an object with xmin, xmax, ymin and ymax")
        try:
            bounds = {key: float(domain[key]) for key in ("xmin", "xmax", "ymin", "ymax")}
        except (KeyError, TypeError, ValueError) as exc:
            raise Stl2DError("domain requires numeric xmin, xmax, ymin and ymax") from exc
        if not all(math.isfinite(value) for value in bounds.values()):
            raise Stl2DError("domain bounds must be finite")
        if bounds["xmin"] >= bounds["xmax"] or bounds["ymin"] >= bounds["ymax"]:
            raise Stl2DError("domain bounds must have positive width and height")
        if not (
            min(p[0] for p in loop) > bounds["xmin"] + tol
            and max(p[0] for p in loop) < bounds["xmax"] - tol
            and min(p[1] for p in loop) > bounds["ymin"] + tol
            and max(p[1] for p in loop) < bounds["ymax"] - tol
        ):
            raise Stl2DError("outer domain must strictly enclose the STL body with the requested tolerance")
        outer = domain.get("boundaries", {})
        required = {"left", "right", "bottom", "top"}
        if not isinstance(outer, dict) or set(outer) != required:
            raise Stl2DError("domain.boundaries must define left, right, bottom and top")
        names = [str(outer[side]).strip() for side in ("left", "right", "bottom", "top")]
        if any(not name for name in names):
            raise Stl2DError("the four outer domain boundaries require non-empty names")
        body_name = str(config.get("body_boundary_name", "Wall")).strip()
        if not body_name or body_name in names:
            raise Stl2DError("body_boundary_name must be non-empty and distinct from outer boundary names")
        config["_weld_tolerance"] = tol
        config["_domain_bounds"] = bounds
        config["_outer_boundaries"] = {side: str(outer[side]).strip() for side in required}
        config["_body_boundary_name"] = body_name
        return config, {}
    entries = config.get("boundaries")
    if not isinstance(entries, list) or not entries:
        raise Stl2DError("sidecar requires boundaries or a domain configuration")
    n, assignments = len(loop), {}
    for entry in entries:
        name, ends = str(entry.get("name", "")).strip(), entry.get("endpoints")
        if not name or not isinstance(ends, list) or len(ends) != 2:
            raise Stl2DError("each boundary requires a name and two [x, y] endpoints")
        indices = []
        for xy in ends:
            if not isinstance(xy, list) or len(xy) != 2:
                raise Stl2DError(f"{name}: endpoints must be [x, y]")
            matches = [i for i,p in enumerate(loop) if math.hypot(p[0]-float(xy[0]),p[1]-float(xy[1])) <= tol]
            if len(matches) != 1:
                raise Stl2DError(f"{name}: endpoint does not uniquely match an STL boundary vertex")
            indices.append(matches[0])
        a,b = indices
        if (a-b) % n not in (1,n-1):
            raise Stl2DError(f"{name}: endpoints must identify one consecutive edge")
        edge = tuple(sorted((a,b)))
        if edge in assignments:
            raise Stl2DError(f"boundary edge assigned more than once: {name}")
        assignments[edge] = name
    expected = {tuple(sorted((i,(i+1)%n))) for i in range(n)}
    if set(assignments) != expected:
        raise Stl2DError(f"assign every outer edge exactly once (missing={len(expected-set(assignments))})")
    config["_weld_tolerance"] = tol
    return config, assignments

def _transfinite_cell_layout(loop: list[Point], config: dict) -> dict | None:
    raw = config.get("transfinite_cells")
    if raw is None and not config.get("one_quad_cell", False):
        return None
    if "_domain_bounds" in config:
        raise Stl2DError("transfinite cell layout is currently supported for a single STL boundary only")
    if raw is None:
        nx, ny = 1, 1
    else:
        if not isinstance(raw, dict):
            raise Stl2DError("transfinite_cells must define integer x and y cell counts")
        nx, ny = raw.get("x"), raw.get("y")
        if any(isinstance(value, bool) or not isinstance(value, int) or value < 1 for value in (nx, ny)):
            raise Stl2DError("transfinite_cells.x and .y must be positive integers")
    if len(loop) != 4:
        raise Stl2DError("transfinite quadrilateral layout requires exactly four STL boundary vertices")
    tolerance = max(float(config.get("_weld_tolerance", 1e-9)), 1e-12)
    horizontal, vertical = [], []
    for index, point in enumerate(loop):
        following = loop[(index + 1) % len(loop)]
        dx, dy = abs(following[0] - point[0]), abs(following[1] - point[1])
        curve = index + 1
        if dy <= tolerance and dx > tolerance:
            horizontal.append(curve)
        elif dx <= tolerance and dy > tolerance:
            vertical.append(curve)
        else:
            raise Stl2DError("transfinite quadrilateral layout requires axis-aligned horizontal and vertical edges")
    if len(horizontal) != 2 or len(vertical) != 2:
        raise Stl2DError("transfinite quadrilateral requires two horizontal and two vertical edges")
    return {"x": nx, "y": ny, "horizontal_curves": horizontal, "vertical_curves": vertical}

def build_geo(loop: list[Point], assignments: dict[tuple[int,int],str], config: dict) -> str:
    size = float(config.get("mesh_size", 1.0))
    if not math.isfinite(size) or size <= 0:
        raise Stl2DError("mesh_size must be positive and finite")
    outer_size = float(config.get("outer_mesh_size", size))
    if "_domain_bounds" in config and (
        not math.isfinite(outer_size) or outer_size <= 0
    ):
        raise Stl2DError("outer_mesh_size must be positive and finite")
    cell_layout = _transfinite_cell_layout(loop, config)
    if cell_layout is None and "_domain_bounds" not in config:
        recombine_surface = config.get("recombine_surface", False)
        if not isinstance(recombine_surface, bool):
            raise Stl2DError("recombine_surface must be a boolean")
        if recombine_surface:
            raise Stl2DError(
                "Automatic triangle meshing does not support recombine_surface"
            )
    lines = ['// Generated from planar STL by OneFLOW.', 'SetFactory("Built-in");', f"lc = {size:.17g};"]
    if cell_layout is None and "_domain_bounds" not in config:
        lines.extend([
            f"Mesh.MeshSizeMin = {size:.17g};",
            f"Mesh.MeshSizeMax = {size:.17g};",
            "Mesh.MeshSizeFromPoints = 1;",
            "Mesh.MeshSizeExtendFromBoundary = 1;",
        ])
    if "_domain_bounds" in config:
        lines.extend([
            f"lc_outer = {outer_size:.17g};",
            f"Mesh.MeshSizeMin = {min(size, outer_size):.17g};",
            f"Mesh.MeshSizeMax = {max(size, outer_size):.17g};",
        ])
        n = len(loop)
        bounds = config["_domain_bounds"]
        body_lines = list(range(1, n + 1))
        for i,(x,y,_z) in enumerate(loop,1):
            lines.append(f"Point({i}) = {{{x:.17g}, {y:.17g}, 0, lc}};")
        # Orient the internal wall curves clockwise, matching the fluid-side
        # boundary orientation; the outer rectangle remains counter-clockwise.
        reverse_loop = [1, *range(n, 1, -1), 1]
        for i in range(n):
            lines.append(f"Line({i+1}) = {{{reverse_loop[i]}, {reverse_loop[i+1]}}};")
        outer_points = [
            (bounds["xmin"], bounds["ymin"]),
            (bounds["xmax"], bounds["ymin"]),
            (bounds["xmax"], bounds["ymax"]),
            (bounds["xmin"], bounds["ymax"]),
        ]
        first_outer_point = n + 1
        for offset,(x,y) in enumerate(outer_points):
            lines.append(f"Point({first_outer_point+offset}) = {{{x:.17g}, {y:.17g}, 0, lc_outer}};")
        first_outer_line = n + 1
        outer_line_ids = {
            "bottom": first_outer_line,
            "right": first_outer_line + 1,
            "top": first_outer_line + 2,
            "left": first_outer_line + 3,
        }
        for offset in range(4):
            a = first_outer_point + offset
            b = first_outer_point + (offset + 1) % 4
            lines.append(f"Line({first_outer_line+offset}) = {{{a}, {b}}};")
        lines.append(f"Curve Loop(1) = {{{', '.join(str(i) for i in outer_line_ids.values())}}};")
        lines.append(f"Curve Loop(2) = {{{', '.join(map(str, body_lines))}}};")
        lines.append("Plane Surface(1) = {1, 2};")
        body_name = re.sub(r"[^A-Za-z0-9_.-]+", "_", config["_body_boundary_name"])
        lines.append(f'Physical Curve("{body_name}") = {{{", ".join(map(str, body_lines))}}};')
        outer_groups = defaultdict(list)
        for side in ("left", "right", "bottom", "top"):
            name = re.sub(r"[^A-Za-z0-9_.-]+", "_", config["_outer_boundaries"][side])
            outer_groups[name].append(outer_line_ids[side])
        for name, curves in sorted(outer_groups.items()):
            lines.append(f'Physical Curve("{name}") = {{{", ".join(map(str, curves))}}};')
    else:
        for i,(x,y,_z) in enumerate(loop,1):
            lines.append(f"Point({i}) = {{{x:.17g}, {y:.17g}, 0, lc}};")
        groups = defaultdict(list)
        for i in range(len(loop)):
            curve = i+1
            lines.append(f"Line({curve}) = {{{i+1}, {(i+1)%len(loop)+1}}};")
            groups[assignments[tuple(sorted((i,(i+1)%len(loop))))]].append(curve)
        lines.extend([f"Curve Loop(1) = {{{', '.join(str(i) for i in range(1,len(loop)+1))}}};", "Plane Surface(1) = {1};"])
        for name,curves in sorted(groups.items()):
            safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", name)
            lines.append(f'Physical Curve("{safe}") = {{{", ".join(map(str,curves))}}};')
        if cell_layout:
            horizontal = ", ".join(map(str, cell_layout["horizontal_curves"]))
            vertical = ", ".join(map(str, cell_layout["vertical_curves"]))
            lines.append(f"Transfinite Curve {{{horizontal}}} = {cell_layout['x'] + 1};")
            lines.append(f"Transfinite Curve {{{vertical}}} = {cell_layout['y'] + 1};")
            lines.append("Transfinite Surface {1};")
            recombine_surface = config.get("recombine_surface", True)
            if not isinstance(recombine_surface, bool):
                raise Stl2DError("recombine_surface must be a boolean")
            if recombine_surface:
                lines.append("Recombine Surface {1};")
    surface = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(config.get("surface_name","Fluid")))
    lines.append(f'Physical Surface("{surface}") = {{1}};')
    lines.extend(["Mesh 2;", ""])
    return "\n".join(lines)

def _cube_volume_geo(
    triangles: list[Triangle],
    tolerance: float,
    boundary_groups: dict[str, str],
    volume_mesh: str,
    mesh_size: float = 0.5,
) -> tuple[str, dict]:
    surface_by_side = {
        "xmin": 3, "xmax": 4, "ymin": 5,
        "ymax": 6, "zmin": 1, "zmax": 2,
    }
    if not isinstance(boundary_groups, dict) or set(boundary_groups) != set(surface_by_side):
        raise Stl2DError(
            "cube volume boundary_groups must define xmin, xmax, ymin, ymax, zmin and zmax"
        )
    physical_groups: dict[str, int] = {}
    for side, surface_id in surface_by_side.items():
        semantic = boundary_groups[side]
        if not isinstance(semantic, str) or semantic.lower() not in {
            "inflow", "outflow", "wall", "farfield", "symmetry"
        }:
            raise Stl2DError(f"cube volume has an unsupported boundary semantic for {side}")
        physical_name = f"{side}_{semantic.lower()}"
        if not re.fullmatch(r"[A-Za-z0-9_.-]{1,32}", physical_name):
            raise Stl2DError(f"cube volume has an invalid Physical Group name: {physical_name}")
        if physical_name in physical_groups:
            raise Stl2DError(f"cube volume Physical Group names must be unique: {physical_name}")
        physical_groups[physical_name] = surface_id

    if len(triangles) != 12:
        raise Stl2DError("cube volume requires exactly 12 STL surface triangles")
    points = [point for triangle in triangles for point in triangle]
    lower = tuple(min(point[axis] for point in points) for axis in range(3))
    upper = tuple(max(point[axis] for point in points) for axis in range(3))
    if any(not math.isfinite(lower[a]) or upper[a] - lower[a] <= tolerance for a in range(3)):
        raise Stl2DError("cube STL must have finite positive extent on all three axes")
    extents = tuple(upper[a] - lower[a] for a in range(3))
    if max(extents) - min(extents) > tolerance:
        raise Stl2DError("cube volume requires equal axis-aligned edge lengths")

    corners: tuple[Point, ...] = (
        (lower[0],lower[1],lower[2]), (upper[0],lower[1],lower[2]),
        (upper[0],upper[1],lower[2]), (lower[0],upper[1],lower[2]),
        (lower[0],lower[1],upper[2]), (upper[0],lower[1],upper[2]),
        (upper[0],upper[1],upper[2]), (lower[0],upper[1],upper[2]))
    def corner_id(point: Point) -> int:
        matches = [i for i, corner in enumerate(corners)
                   if all(abs(point[a] - corner[a]) <= tolerance for a in range(3))]
        if len(matches) != 1:
            raise Stl2DError("STL vertex is not a unique axis-aligned cube corner")
        return matches[0]
    face_counts = [0] * 6
    edge_counts: Counter[tuple[int, int]] = Counter()
    for triangle in triangles:
        ids = [corner_id(point) for point in triangle]
        face = None
        for axis in range(3):
            value = triangle[0][axis]
            if all(abs(point[axis] - value) <= tolerance for point in triangle):
                if abs(value - lower[axis]) <= tolerance: face = 2 * axis
                elif abs(value - upper[axis]) <= tolerance: face = 2 * axis + 1
        if face is None:
            raise Stl2DError("each STL triangle must lie on one cube bounding-box face")
        face_counts[face] += 1
        for a, b in ((ids[0],ids[1]),(ids[1],ids[2]),(ids[2],ids[0])):
            edge_counts[tuple(sorted((a,b)))] += 1
    if face_counts != [2] * 6:
        raise Stl2DError("each cube face must contain exactly two STL triangles")
    if len(edge_counts) != 18 or any(count != 2 for count in edge_counts.values()):
        raise Stl2DError("cube STL must form a closed two-manifold surface")
    lines = ['SetFactory("Built-in");', "Mesh.MshFileVersion = 4.1;",
             "Mesh.Binary = 0;", "Mesh.ElementOrder = 1;"]
    for i, point in enumerate(corners, 1):
        lines.append(f"Point({i}) = {{{point[0]:.17g}, {point[1]:.17g}, {point[2]:.17g}, 1}};")
    edges = ((1,2),(2,3),(3,4),(4,1),(5,6),(6,7),(7,8),(8,5),(1,5),(2,6),(3,7),(4,8))
    for i, (a,b) in enumerate(edges, 1):
        lines.append(f"Line({i}) = {{{a}, {b}}};")
    loops = ("{1, 2, 3, 4}", "{5, 6, 7, 8}", "{4, 9, -8, -12}",
             "{2, 11, -6, -10}", "{1, 10, -5, -9}", "{3, 12, -7, -11}")
    for i, loop in enumerate(loops, 1):
        lines.extend((f"Curve Loop({i}) = {loop};", f"Plane Surface({i}) = {{{i}}};"))
    lines.extend(("Surface Loop(1) = {-1, 2, 3, 4, 5, 6};", "Volume(1) = {1};"))
    if volume_mesh == "axis_aligned_cube_single_hex":
        lines.extend(("Transfinite Curve {1:12} = 2;", "Transfinite Surface {1:6};",
                      "Recombine Surface {1:6};", "Transfinite Volume {1};",
                      "Recombine Volume {1};"))
    elif volume_mesh == "axis_aligned_cube_six_tetrahedra":
        lines.extend((
            "Transfinite Curve {1:12} = 2;",
            "Transfinite Surface {1} = {1, 2, 3, 4};",
            "Transfinite Surface {2} = {5, 6, 7, 8};",
            "Transfinite Surface {3} = {4, 1, 5, 8} Right;",
            "Transfinite Surface {4} = {2, 3, 7, 6};",
            "Transfinite Surface {5} = {1, 2, 6, 5};",
            "Transfinite Surface {6} = {3, 4, 8, 7} Right;",
            "Transfinite Volume {1} = {1, 2, 3, 4, 5, 6, 7, 8};",
        ))
    elif volume_mesh == "axis_aligned_cube_automatic_tetrahedra":
        if not math.isfinite(mesh_size) or mesh_size <= 0:
            raise Stl2DError("mesh_size must be positive and finite")
        lines.extend((f"Mesh.MeshSizeMin = {mesh_size:.17g};",
                      f"Mesh.MeshSizeMax = {mesh_size:.17g};",
                      "Mesh.MeshSizeFromPoints = 1;",
                      "Mesh.MeshSizeExtendFromBoundary = 1;"))
    else:
        raise Stl2DError(f"unsupported cube volume mesh strategy: {volume_mesh}")
    for name, surface_id in physical_groups.items():
        lines.append(f'Physical Surface("{name}") = {{{surface_id}}};')
    lines.append('Physical Volume("fluid_domain") = {1};')
    return "\n".join(lines) + "\n", {
        "triangles": 12, "unique_cube_corners": 8,
        "boundary_faces": 6, "volume_cells": 1,
        "physical_groups": sorted([*physical_groups, "fluid_domain"]),
    }


def convert(stl: Path, sidecar: Path, geo: Path, tolerance: float = 1e-9) -> dict:
    triangles = read_stl(stl)
    sidecar_data = json.loads(sidecar.read_text(encoding="utf-8"))
    volume_mesh = sidecar_data.get("volume_mesh")
    supported_volume_meshes = {
        "axis_aligned_cube_single_hex",
        "axis_aligned_cube_six_tetrahedra",
        "axis_aligned_cube_automatic_tetrahedra",
    }
    if sidecar_data.get("mesh_dimension") == 3 and volume_mesh not in supported_volume_meshes:
        raise Stl2DError(f"unsupported 3D cube volume mesh strategy: {volume_mesh}")
    if volume_mesh in {
        "axis_aligned_cube_single_hex",
        "axis_aligned_cube_six_tetrahedra",
        "axis_aligned_cube_automatic_tetrahedra",
    }:
        geometry, report = _cube_volume_geo(
            triangles, tolerance, sidecar_data.get("boundary_groups"),
            volume_mesh, float(sidecar_data.get("mesh_size", 0.5)))
        geo.parent.mkdir(parents=True, exist_ok=True)
        geo.write_text(geometry, encoding="utf-8", newline="\n")
        report["geo"] = str(geo)
        return report
    loop = extract_boundary(triangles, tolerance)
    config, assignments = _read_sidecar(sidecar, loop, tolerance)
    layout = _transfinite_cell_layout(loop, config)
    geo.parent.mkdir(parents=True, exist_ok=True)
    geo.write_text(build_geo(loop, assignments, config), encoding="utf-8", newline="\n")
    report = {
        "triangles": len(triangles),
        "boundary_vertices": len(loop),
        "boundary_edges": len(assignments),
        "body_boundary_edges": len(loop) if "_domain_bounds" in config else 0,
        "geo": str(geo),
    }
    if layout:
        report["structured_cells"] = {
            "x": layout["x"],
            "y": layout["y"],
            "total": layout["x"] * layout["y"],
        }
    return report


def _validate_tetrahedral_boundary(elements: list[tuple[int, int, tuple[int, ...]]]) -> dict:
    """Check TRI_3 facets against the faces of an arbitrary linear tetrahedral mesh."""
    face_incidence: Counter[tuple[int, int, int]] = Counter()
    surface_faces: set[tuple[int, int, int]] = set()
    tetrahedron_count = 0
    for entity_dimension, element_type, node_ids in elements:
        if entity_dimension == 3:
            if element_type != 4 or len(node_ids) != 4:
                raise Stl2DError(
                    "tetrahedral topology validation requires first-order TETRA_4 volume cells"
                )
            if len(set(node_ids)) != 4:
                raise Stl2DError("TETRA_4 element contains repeated node tags")
            tetrahedron_count += 1
            for omitted in range(4):
                face = tuple(sorted(node_ids[:omitted] + node_ids[omitted + 1:]))
                face_incidence[face] += 1
                if face_incidence[face] > 2:
                    raise Stl2DError(
                        f"non-manifold tetrahedral face is shared by more than two cells: {face}"
                    )
        elif entity_dimension == 2:
            if element_type != 2 or len(node_ids) != 3:
                raise Stl2DError(
                    "tetrahedral boundary validation requires first-order TRI_3 surface elements"
                )
            face = tuple(sorted(node_ids))
            if face in surface_faces:
                raise Stl2DError(f"duplicate TRI_3 boundary element: {face}")
            surface_faces.add(face)

    if tetrahedron_count == 0:
        raise Stl2DError("tetrahedral topology validation found no TETRA_4 volume cells")
    for face in surface_faces:
        if face not in face_incidence:
            raise Stl2DError(
                f"TRI_3 surface element does not match any TETRA_4 face: {face}"
            )
    boundary_faces = {face for face, count in face_incidence.items() if count == 1}
    missing = boundary_faces - surface_faces
    if missing:
        face = min(missing)
        raise Stl2DError(
            f"TETRA_4 exterior face has no matching TRI_3 boundary element: {face}"
        )
    return {
        "tetrahedron_cells": tetrahedron_count,
        "exterior_faces": len(boundary_faces),
        "internal_faces": sum(count == 2 for count in face_incidence.values()),
        "surface_triangles": len(surface_faces),
    }


def inspect_msh4(path: Path) -> dict:
    """Inspect MSH 4.1 topology, coordinates and first-order triangle quality."""
    lines = path.read_text(encoding="utf-8").splitlines()

    def section(name: str) -> list[str]:
        try:
            first = lines.index("$" + name)
            last = lines.index("$End" + name, first + 1)
        except ValueError as exc:
            raise Stl2DError(f"MSH file is missing the {name} section") from exc
        return lines[first + 1:last]

    mesh_format = section("MeshFormat")
    if not mesh_format:
        raise Stl2DError("MSH MeshFormat section is empty")
    format_values = mesh_format[0].split()
    if len(format_values) < 3 or format_values[0] != "4.1" or format_values[1] != "0":
        raise Stl2DError("generated mesh must be MSH 4.1 ASCII")

    node_tokens = " ".join(section("Nodes")).split()
    if len(node_tokens) < 4:
        raise Stl2DError("MSH Nodes section is empty or truncated")
    block_count, node_count = int(node_tokens[0]), int(node_tokens[1])
    cursor = 4
    nodes: dict[int, Point] = {}
    for _ in range(block_count):
        if cursor + 4 > len(node_tokens):
            raise Stl2DError("truncated MSH 4.1 node block header")
        entity_dim, _entity_tag, parametric, count = map(int, node_tokens[cursor:cursor + 4])
        cursor += 4
        if entity_dim not in (0, 1, 2, 3) or parametric not in (0, 1) or count < 0:
            raise Stl2DError("invalid MSH 4.1 node block header")
        tags = list(map(int, node_tokens[cursor:cursor + count]))
        if len(tags) != count:
            raise Stl2DError("truncated MSH 4.1 node tags")
        cursor += count
        for tag in tags:
            coordinate_count = 3 + (entity_dim if parametric else 0)
            if cursor + coordinate_count > len(node_tokens):
                raise Stl2DError("truncated MSH 4.1 node coordinates")
            xyz = tuple(float(value) for value in node_tokens[cursor:cursor + 3])
            cursor += coordinate_count
            if not all(math.isfinite(value) for value in xyz) or tag in nodes:
                raise Stl2DError("non-finite coordinates or duplicate MSH node tag")
            nodes[tag] = xyz
    if len(nodes) != node_count:
        raise Stl2DError("MSH node block counts do not match the file header")

    coordinate_values = list(nodes.values())
    coordinate_scale = max(
        max(point[axis] for point in coordinate_values)
        - min(point[axis] for point in coordinate_values)
        for axis in range(3)
    )
    area_tolerance = 64.0 * math.ulp(1.0) * coordinate_scale * coordinate_scale

    element_tokens = " ".join(section("Elements")).split()
    if len(element_tokens) < 4:
        raise Stl2DError("MSH Elements section is empty or truncated")
    element_block_count, element_count = int(element_tokens[0]), int(element_tokens[1])
    cursor = 4
    element_types: Counter[int] = Counter()
    element_records: list[tuple[int, int, tuple[int, ...]]] = []
    triangle_areas: list[float] = []
    triangle_qualities: list[float] = []
    triangle_min_angles: list[float] = []
    positive_xy = negative_xy = 0
    for _ in range(element_block_count):
        if cursor + 4 > len(element_tokens):
            raise Stl2DError("truncated MSH 4.1 element block header")
        entity_dim, _entity_tag, element_type, count = map(
            int, element_tokens[cursor:cursor + 4]
        )
        cursor += 4
        arity = {1: 2, 2: 3, 3: 4, 4: 4, 5: 8}.get(element_type)
        if arity is None:
            raise Stl2DError(f"unsupported element type {element_type} in reference mesh")
        element_types[element_type] += count
        for _element in range(count):
            if cursor + 1 + arity > len(element_tokens):
                raise Stl2DError("truncated MSH 4.1 element connectivity")
            _tag = int(element_tokens[cursor])
            node_ids = list(map(int, element_tokens[cursor + 1:cursor + 1 + arity]))
            cursor += 1 + arity
            if any(node_id not in nodes for node_id in node_ids):
                raise Stl2DError("MSH element references an unknown node")
            element_records.append((entity_dim, element_type, tuple(node_ids)))
            if element_type != 2:
                continue
            p0, p1, p2 = (nodes[node_id] for node_id in node_ids)
            u = tuple(p1[i] - p0[i] for i in range(3))
            v = tuple(p2[i] - p0[i] for i in range(3))
            cross = (
                u[1] * v[2] - u[2] * v[1],
                u[2] * v[0] - u[0] * v[2],
                u[0] * v[1] - u[1] * v[0],
            )
            area = 0.5 * math.sqrt(sum(value * value for value in cross))
            signed_xy = 0.5 * (u[0] * v[1] - u[1] * v[0])
            lengths = (
                math.dist(p0, p1), math.dist(p1, p2), math.dist(p2, p0)
            )
            if area <= area_tolerance or min(lengths) <= 0.0:
                raise Stl2DError("generated mesh contains a degenerate TRI_3 element")
            angles = []
            for opposite, adjacent_a, adjacent_b in (
                (lengths[1], lengths[0], lengths[2]),
                (lengths[2], lengths[0], lengths[1]),
                (lengths[0], lengths[1], lengths[2]),
            ):
                cosine = (adjacent_a**2 + adjacent_b**2 - opposite**2) / (
                    2.0 * adjacent_a * adjacent_b
                )
                angles.append(math.degrees(math.acos(max(-1.0, min(1.0, cosine)))))
            triangle_areas.append(area)
            triangle_qualities.append(
                4.0 * math.sqrt(3.0) * area / sum(length * length for length in lengths)
            )
            triangle_min_angles.append(min(angles))
            if signed_xy > 0.0:
                positive_xy += 1
            elif signed_xy < 0.0:
                negative_xy += 1

    if sum(element_types.values()) != element_count:
        raise Stl2DError("MSH 4.1 element block counts do not match the file header")

    tetrahedral_topology = None
    if element_types.get(4, 0):
        volume_types = {
            element_type for dimension, element_type, _nodes in element_records
            if dimension == 3
        }
        if volume_types != {4}:
            raise Stl2DError(
                "automatic tetrahedral boundary validation does not accept mixed volume cell families"
            )
        tetrahedral_topology = _validate_tetrahedral_boundary(element_records)

    physical_groups = []
    physical_section = section("PhysicalNames")
    if physical_section:
        expected_groups = int(physical_section[0])
        for row in physical_section[1:]:
            values = row.split(maxsplit=2)
            if len(values) != 3:
                raise Stl2DError("invalid MSH PhysicalNames record")
            physical_groups.append({
                "dimension": int(values[0]),
                "tag": int(values[1]),
                "name": values[2].strip().strip('"'),
            })
        if len(physical_groups) != expected_groups:
            raise Stl2DError("MSH PhysicalNames count does not match its header")

    coordinates = list(nodes.values())
    bounds = {
        "minimum": [min(point[axis] for point in coordinates) for axis in range(3)],
        "maximum": [max(point[axis] for point in coordinates) for axis in range(3)],
    }
    return {
        "format": format_values[0],
        "binary": False,
        "nodes": node_count,
        "elements": element_count,
        "element_types": {str(key): value for key, value in sorted(element_types.items())},
        "line_elements": element_types.get(1, 0),
        "triangle_elements": element_types.get(2, 0),
        "quadrilateral_elements": element_types.get(3, 0),
        "tetrahedron_elements": element_types.get(4, 0),
        "hexahedron_elements": element_types.get(5, 0),
        "tetrahedral_topology": tetrahedral_topology,
        "physical_groups": physical_groups,
        "coordinate_bounds": bounds,
        "triangle_quality": {
            "degenerate_count": 0,
            "positive_xy_orientation_count": positive_xy,
            "negative_xy_orientation_count": negative_xy,
            "minimum_area": min(triangle_areas) if triangle_areas else 0.0,
            "maximum_area": max(triangle_areas) if triangle_areas else 0.0,
            "minimum_quality": min(triangle_qualities) if triangle_qualities else 0.0,
            "minimum_angle_degrees": min(triangle_min_angles) if triangle_min_angles else 0.0,
        },
    }

def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stl", type=Path)
    parser.add_argument("--boundaries", type=Path, required=True)
    parser.add_argument("--geo", type=Path, required=True)
    parser.add_argument("--tolerance", type=float, default=1e-9)
    parser.add_argument("--gmsh", help="optional Gmsh executable")
    parser.add_argument("--msh", type=Path, help="output MSH path; required with --gmsh")
    args = parser.parse_args()
    try:
        report = convert(args.stl,args.boundaries,args.geo,args.tolerance)
        print(json.dumps(report,indent=2))
        if args.gmsh:
            if not args.msh: parser.error("--msh is required with --gmsh")
            args.msh.parent.mkdir(parents=True, exist_ok=True)
            subprocess.run(
                [args.gmsh, str(args.geo), "-2", "-format", "msh4", "-o", str(args.msh)],
                check=True,
            )
    except (OSError,json.JSONDecodeError,Stl2DError,subprocess.CalledProcessError) as exc:
        parser.exit(2,f"stl_to_gmsh: {exc}\n")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
