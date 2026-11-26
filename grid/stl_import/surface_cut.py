"""Generic plane cutting for closed STL triangle surfaces.

This module deliberately knows nothing about cylinders, radii, centerlines or
analytic surface equations. It clips an input triangle shell against two
user-supplied planes, extracts the resulting boundary loops and triangulates
the planar end caps.
"""
from __future__ import annotations

import math
from collections import Counter, defaultdict
from typing import Iterable

from grid.stl_import.reader import Point, Triangle
from grid.stl_import.surface import StlGeometryError, _weld_vertices


def _sub(a: Point, b: Point) -> Point:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _dot(a: Point, b: Point) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _cross(a: Point, b: Point) -> Point:
    return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])


def _norm(vector: Point) -> float:
    return math.sqrt(_dot(vector, vector))


def _unit(vector: Point, label: str) -> Point:
    length = _norm(vector)
    if not math.isfinite(length) or length <= 0.0:
        raise StlGeometryError(f"{label} must be finite and non-zero")
    return (vector[0]/length, vector[1]/length, vector[2]/length)


def _plane_distance(point: Point, origin: Point, normal: Point) -> float:
    return _dot(_sub(point, origin), normal)


def _project(point: Point, origin: Point, normal: Point) -> Point:
    distance = _plane_distance(point, origin, normal)
    return (point[0]-distance*normal[0], point[1]-distance*normal[1],
            point[2]-distance*normal[2])


def _clean(points: list[Point], tolerance: float) -> list[Point]:
    result = []
    for point in points:
        if not result or math.dist(point, result[-1]) > tolerance:
            result.append(point)
    if len(result) > 1 and math.dist(result[0], result[-1]) <= tolerance:
        result.pop()
    return result


def _clip_polygon(polygon: list[Point], origin: Point, normal: Point,
                  keep_positive: bool, tolerance: float) -> list[Point]:
    if not polygon:
        return []
    def signed(point):
        return _plane_distance(point, origin, normal)
    def inside(distance):
        return distance >= -tolerance if keep_positive else distance <= tolerance

    result = []
    previous = polygon[-1]
    previous_distance = signed(previous)
    previous_inside = inside(previous_distance)
    for current in polygon:
        current_distance = signed(current)
        current_inside = inside(current_distance)
        if current_inside != previous_inside:
            denominator = previous_distance - current_distance
            if abs(denominator) <= 1.0e-30:
                raise StlGeometryError("unstable triangle-plane intersection")
            fraction = previous_distance / denominator
            intersection = tuple(
                previous[axis] + fraction*(current[axis]-previous[axis])
                for axis in range(3)
            )
            result.append(_project(intersection, origin, normal))
        if current_inside:
            result.append(
                _project(current, origin, normal)
                if abs(current_distance) <= tolerance else current
            )
        previous, previous_distance, previous_inside = (
            current, current_distance, current_inside
        )
    return _clean(result, tolerance)


def _validate_shell(triangles: list[Triangle], tolerance: float) -> None:
    if not triangles:
        raise StlGeometryError("STL contains no triangles")
    raw = [point for triangle in triangles for point in triangle]
    points, triangle_ids = _weld_vertices(raw, tolerance)
    uses = defaultdict(list)
    adjacency = defaultdict(set)
    signed_volume = 0.0
    for triangle, ids in zip(triangles, triangle_ids):
        if len(set(ids)) != 3:
            raise StlGeometryError("STL contains a degenerate triangle after vertex welding")
        a, b, c = triangle
        signed_volume += (
            a[0] * (b[1]*c[2] - b[2]*c[1])
            + a[1] * (b[2]*c[0] - b[0]*c[2])
            + a[2] * (b[0]*c[1] - b[1]*c[0])
        ) / 6.0
        for index in range(3):
            first, second = ids[index], ids[(index+1) % 3]
            uses[tuple(sorted((first, second)))].append(
                1 if first < second else -1
            )
            adjacency[first].add(second)
            adjacency[second].add(first)
    if any(len(values) != 2 or sum(values) != 0 for values in uses.values()):
        raise StlGeometryError("generic cutting requires one closed, consistently oriented, two-manifold STL shell")
    if not math.isfinite(signed_volume) or signed_volume <= tolerance**3:
        raise StlGeometryError("STL triangle winding must point outward and enclose positive volume")
    pending, visited = [0], set()
    while pending:
        current = pending.pop()
        if current in visited:
            continue
        visited.add(current)
        pending.extend(adjacency[current] - visited)
    if len(visited) != len(points):
        raise StlGeometryError("generic cutting requires one connected STL shell")


def _extract_loop(edges: list[tuple[int, int]], points: list[Point], label: str) -> list[Point]:
    if len(edges) < 3:
        raise StlGeometryError(f"{label} cut did not produce a closed boundary loop")
    adjacency = defaultdict(set)
    for a, b in edges:
        adjacency[a].add(b)
        adjacency[b].add(a)
    if any(len(neighbors) != 2 for neighbors in adjacency.values()):
        raise StlGeometryError(f"{label} cut boundary is open, branched, or self-touching")
    start, current, previous = min(adjacency), min(adjacency), -1
    ordered = [start]
    while True:
        candidates = adjacency[current] - ({previous} if previous >= 0 else set())
        following = min(candidates)
        if following == start:
            break
        if following in ordered:
            raise StlGeometryError(f"{label} cut boundary repeats a vertex")
        ordered.append(following)
        previous, current = current, following
        if len(ordered) > len(adjacency):
            raise StlGeometryError(f"{label} cut boundary traversal did not close")
    if len(ordered) != len(adjacency):
        raise StlGeometryError(f"{label} cut produced multiple boundary loops")
    return [points[index] for index in ordered]


def _orient2(a, b, c):
    return (b[0]-a[0])*(c[1]-a[1]) - (b[1]-a[1])*(c[0]-a[0])


def _triangulate_cap(loop: list[Point], plane_normal: Point,
                     outward_normal: Point, tolerance: float) -> list[Triangle]:
    reference = min(((1.,0.,0.), (0.,1.,0.), (0.,0.,1.)),
                    key=lambda axis: abs(_dot(axis, plane_normal)))
    u = _unit(_cross(reference, plane_normal), "projection basis")
    v = _cross(plane_normal, u)
    origin = loop[0]
    points2 = [(_dot(_sub(p, origin), u), _dot(_sub(p, origin), v)) for p in loop]
    area2 = sum(
        points2[i][0]*points2[(i+1)%len(points2)][1]
        - points2[(i+1)%len(points2)][0]*points2[i][1]
        for i in range(len(points2))
    )
    if abs(area2) <= tolerance*tolerance:
        raise StlGeometryError("cut boundary encloses zero area")
    winding = 1.0 if area2 > 0 else -1.0
    remaining = list(range(len(loop)))
    triangles = []
    area_tolerance = tolerance*tolerance

    def contains(p, a, b, c):
        return (winding*_orient2(a,b,p) >= -area_tolerance
                and winding*_orient2(b,c,p) >= -area_tolerance
                and winding*_orient2(c,a,p) >= -area_tolerance)

    while len(remaining) > 3:
        found = False
        for position, current in enumerate(remaining):
            previous = remaining[position-1]
            following = remaining[(position+1) % len(remaining)]
            a, b, c = points2[previous], points2[current], points2[following]
            if winding*_orient2(a,b,c) <= area_tolerance:
                continue
            if any(contains(points2[other],a,b,c)
                   for other in remaining if other not in (previous,current,following)):
                continue
            triangles.append((loop[previous],loop[current],loop[following]))
            del remaining[position]
            found = True
            break
        if not found:
            raise StlGeometryError("cut boundary could not be triangulated as a simple planar polygon")
    triangles.append(tuple(loop[index] for index in remaining))
    oriented = []
    for a, b, c in triangles:
        if _dot(_cross(_sub(b,a),_sub(c,a)),outward_normal) < 0:
            b, c = c, b
        oriented.append((a,b,c))
    return oriented


def cut_closed_surface_between_planes(
    triangles: Iterable[Triangle],
    inlet_point: Point,
    inlet_normal: Point,
    outlet_point: Point,
    outlet_normal: Point,
    tolerance: float = 1.0e-8,
    inlet_keep_positive: bool = True,
    outlet_keep_positive: bool = False,
) -> dict:
    """Clip a closed STL shell between two planes and generate its end caps.

    Retained half-spaces are explicit for each plane. Exactly one simple loop
    per plane is supported. No geometry manifest or analytic primitive is read.
    """
    source = list(triangles)
    if not math.isfinite(tolerance) or tolerance <= 0:
        raise StlGeometryError("tolerance must be positive and finite")
    if any(not math.isfinite(value) for point in (inlet_point,outlet_point) for value in point):
        raise StlGeometryError("cut plane points must contain finite coordinates")
    inlet_normal = _unit(inlet_normal, "inlet normal")
    outlet_normal = _unit(outlet_normal, "outlet normal")
    _validate_shell(source, tolerance)

    clipped = []
    for triangle in source:
        polygon = _clip_polygon(
            list(triangle), inlet_point, inlet_normal, inlet_keep_positive, tolerance
        )
        polygon = _clip_polygon(
            polygon, outlet_point, outlet_normal, outlet_keep_positive, tolerance
        )
        polygon = _clean(polygon, tolerance)
        for index in range(1, len(polygon)-1):
            candidate = (polygon[0], polygon[index], polygon[index+1])
            if _norm(_cross(_sub(candidate[1],candidate[0]),
                            _sub(candidate[2],candidate[0]))) > tolerance*tolerance:
                clipped.append(candidate)
    if not clipped:
        raise StlGeometryError("cutting planes removed the complete STL surface")

    raw = [point for triangle in clipped for point in triangle]
    points, triangle_ids = _weld_vertices(raw, tolerance)
    edge_counts = Counter()
    for ids in triangle_ids:
        for index in range(3):
            edge_counts[tuple(sorted((ids[index],ids[(index+1)%3])))] += 1
    inlet_edges, outlet_edges = [], []
    for edge, count in edge_counts.items():
        if count == 2:
            continue
        if count != 1:
            raise StlGeometryError("clipped STL contains a non-manifold edge")
        a, b = (points[index] for index in edge)
        on_inlet = all(abs(_plane_distance(p,inlet_point,inlet_normal)) <= 2*tolerance for p in (a,b))
        on_outlet = all(abs(_plane_distance(p,outlet_point,outlet_normal)) <= 2*tolerance for p in (a,b))
        if on_inlet and not on_outlet:
            inlet_edges.append(edge)
        elif on_outlet and not on_inlet:
            outlet_edges.append(edge)
        else:
            raise StlGeometryError("clipped STL has an unexpected open edge away from a unique cutting plane")

    inlet_loop = _extract_loop(inlet_edges, points, "inlet")
    outlet_loop = _extract_loop(outlet_edges, points, "outlet")
    inlet_outward = (
        tuple(-value for value in inlet_normal)
        if inlet_keep_positive else inlet_normal
    )
    outlet_outward = (
        tuple(-value for value in outlet_normal)
        if outlet_keep_positive else outlet_normal
    )
    inlet_cap = _triangulate_cap(inlet_loop, inlet_normal, inlet_outward, tolerance)
    outlet_cap = _triangulate_cap(outlet_loop, outlet_normal, outlet_outward, tolerance)
    closed = [*clipped, *inlet_cap, *outlet_cap]
    return {
        "side_triangles": clipped,
        "inlet_loop": inlet_loop,
        "outlet_loop": outlet_loop,
        "inlet_cap_triangles": inlet_cap,
        "outlet_cap_triangles": outlet_cap,
        "closed_triangles": closed,
        "source_triangle_count": len(source),
        "side_triangle_count": len(clipped),
        "inlet_loop_vertex_count": len(inlet_loop),
        "outlet_loop_vertex_count": len(outlet_loop),
        "inlet_cap_triangle_count": len(inlet_cap),
        "outlet_cap_triangle_count": len(outlet_cap),
    }
