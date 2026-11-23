"""Non-destructive, scale-aware diagnostics for arbitrary STL surfaces."""
from __future__ import annotations
import math
from collections import defaultdict
from typing import Iterable
from grid.stl_import.reader import Triangle
from grid.stl_import.surface import StlGeometryError, _weld_vertices


def _components(face_count, edge_faces):
    adjacency = [set() for _ in range(face_count)]
    for faces in edge_faces.values():
        for face in faces:
            adjacency[face].update(other for other in faces if other != face)
    pending, result = set(range(face_count)), []
    while pending:
        start = min(pending)
        pending.remove(start)
        stack, component = [start], []
        while stack:
            face = stack.pop()
            component.append(face)
            for other in adjacency[face] & pending:
                pending.remove(other)
                stack.append(other)
        result.append(sorted(component))
    return result


def _boundary_summary(edges):
    adjacency = defaultdict(set)
    for a, b in edges:
        adjacency[a].add(b)
        adjacency[b].add(a)
    pending, component_sizes = set(adjacency), []
    branches = sum(len(neighbors) != 2 for neighbors in adjacency.values())
    while pending:
        start = min(pending)
        pending.remove(start)
        stack, vertex_count = [start], 0
        while stack:
            current = stack.pop()
            vertex_count += 1
            for other in adjacency[current] & pending:
                pending.remove(other)
                stack.append(other)
        component_sizes.append(vertex_count)
    return {
        "boundary_edge_count": len(edges),
        "boundary_vertex_count": len(adjacency),
        "boundary_component_count": len(component_sizes),
        "boundary_branch_vertex_count": branches,
        "boundary_loop_count": len(component_sizes) if edges and branches == 0 else 0,
        "boundary_component_vertex_counts": component_sizes,
        "boundary_is_loop_like": bool(edges) and branches == 0,
    }


def analyze_surface(triangles: Iterable[Triangle], weld_tolerance: float | None = None) -> dict:
    """Inspect a triangle surface without repairing or changing its input."""
    faces = list(triangles)
    if not faces:
        raise StlGeometryError("STL contains no triangles")
    raw = [point for face in faces for point in face]
    if any(len(point) != 3 or not all(math.isfinite(v) for v in point) for point in raw):
        raise StlGeometryError("STL contains a non-finite or invalid coordinate")
    bounds = {
        axis: [min(p[i] for p in raw), max(p[i] for p in raw)]
        for i, axis in enumerate(("x", "y", "z"))
    }
    diagonal = math.sqrt(sum((bounds[a][1] - bounds[a][0]) ** 2 for a in ("x", "y", "z")))
    if not math.isfinite(diagonal):
        raise StlGeometryError("STL coordinate range is too large")
    tolerance = max(diagonal * 1.0e-9, 1.0e-12) if weld_tolerance is None else weld_tolerance
    if not math.isfinite(tolerance) or tolerance <= 0:
        raise StlGeometryError("weld tolerance must be positive and finite")

    points, face_ids = _weld_vertices(raw, tolerance)
    edge_faces, edge_directions = defaultdict(list), defaultdict(list)
    duplicate_faces, degenerate, zero_area, areas, volumes = {}, [], [], [], []
    for fi, (face, ids) in enumerate(zip(faces, face_ids)):
        if len(set(ids)) != 3:
            degenerate.append(fi)
            areas.append(0.0)
            volumes.append(0.0)
            continue
        a, b, c = face
        ab = tuple(b[i] - a[i] for i in range(3))
        ac = tuple(c[i] - a[i] for i in range(3))
        cross = (ab[1]*ac[2]-ab[2]*ac[1], ab[2]*ac[0]-ab[0]*ac[2], ab[0]*ac[1]-ab[1]*ac[0])
        area = 0.5 * math.sqrt(sum(v*v for v in cross))
        areas.append(area)
        if not math.isfinite(area) or area <= tolerance*tolerance:
            zero_area.append(fi)
            degenerate.append(fi)
        key = tuple(sorted(ids))
        duplicate_faces[key] = duplicate_faces.get(key, 0) + 1
        for i in range(3):
            u, v = ids[i], ids[(i+1) % 3]
            edge = tuple(sorted((u, v)))
            edge_faces[edge].append(fi)
            edge_directions[edge].append(1 if u < v else -1)
        volumes.append((a[0]*(b[1]*c[2]-b[2]*c[1])
                      + a[1]*(b[2]*c[0]-b[0]*c[2])
                      + a[2]*(b[0]*c[1]-b[1]*c[0])) / 6.0)

    duplicates = sum(n-1 for n in duplicate_faces.values() if n > 1)
    boundary_edges = [e for e, uses in edge_faces.items() if len(uses) == 1]
    non_manifold = [e for e, uses in edge_faces.items() if len(uses) > 2]
    conflicts = [e for e, uses in edge_faces.items()
                 if len(uses) == 2 and sum(edge_directions[e]) != 0]
    components = _components(len(faces), edge_faces)
    boundary = _boundary_summary(boundary_edges)
    edge_manifold = not non_manifold and not degenerate and not zero_area
    closed = edge_manifold and not boundary_edges
    consistently_oriented = not conflicts
    component_reports = []
    for ci, face_indices in enumerate(components):
        face_set = set(face_indices)
        edges = {e for e, uses in edge_faces.items() if any(fi in face_set for fi in uses)}
        vertices = {v for edge in edges for v in edge}
        is_closed = bool(edges) and all(len(edge_faces[e]) == 2 for e in edges)
        is_oriented = all(len(edge_faces[e]) != 2 or sum(edge_directions[e]) == 0 for e in edges)
        signed_volume = sum(volumes[i] for i in face_indices)
        euler = len(vertices) - len(edges) + len(face_indices)
        genus = ((2-euler)//2 if is_closed and is_oriented and
                 2-euler >= 0 and (2-euler) % 2 == 0 else None)
        outward = signed_volume > tolerance**3 if is_closed and is_oriented else None
        component_reports.append({
            "component": ci, "triangle_count": len(face_indices),
            "vertex_count": len(vertices), "edge_count": len(edges),
            "closed": is_closed, "consistently_oriented": is_oriented,
            "euler_characteristic": euler if is_closed and is_oriented else None,
            "genus": genus, "signed_volume": signed_volume if is_closed and is_oriented else None,
            "outward_winding": outward,
        })

    valid_geometry = not degenerate and not zero_area and duplicates == 0
    if not valid_geometry or non_manifold or conflicts:
        status = "invalid"
    elif boundary_edges or len(components) != 1:
        status = "warning"
    else:
        status = "valid"
    return {
        "status": status, "triangle_count": len(faces),
        "welded_vertex_count": len(points), "unique_edge_count": len(edge_faces),
        "coordinate_bounds": bounds, "bounding_box_diagonal": diagonal,
        "weld_tolerance": tolerance, "surface_area": sum(areas),
        "degenerate_triangle_count": len(degenerate), "zero_area_triangle_count": len(zero_area),
        "duplicate_triangle_count": duplicates, "boundary_edge_count": len(boundary_edges),
        "non_manifold_edge_count": len(non_manifold), "orientation_conflict_edge_count": len(conflicts),
        "connected_component_count": len(components), "closed": closed,
        "edge_manifold": edge_manifold, "consistently_oriented": consistently_oriented,
        "boundary": boundary, "components": component_reports,
        "self_intersection_check": "not_performed", "units": "unknown", "repairs_performed": False,
        "volume_mesh_ready": bool(status == "valid" and closed and consistently_oriented
                                   and len(components) == 1
                                   and component_reports[0]["outward_winding"] is True),
    }
