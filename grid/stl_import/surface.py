"""Topology operations on STL triangle surfaces, independent of mesh generators."""
from __future__ import annotations
import math
from collections import Counter, defaultdict
from typing import Iterable
from grid.stl_import.reader import Point, Triangle

class StlGeometryError(ValueError):
    """Unsupported or invalid STL surface geometry."""

def _check_simple_loop(loop: list[Point], tol: float) -> None:
    def orient(a, b, c):
        return (b[0]-a[0])*(c[1]-a[1]) - (b[1]-a[1])*(c[0]-a[0])
    n = len(loop)
    for i in range(n):
        a, b = loop[i], loop[(i+1) % n]
        for j in range(i+1, n):
            if j == i+1 or (j+1) % n == i:
                continue
            c, d = loop[j], loop[(j+1) % n]
            o1, o2, o3, o4 = orient(a,b,c), orient(a,b,d), orient(c,d,a), orient(c,d,b)
            if ((o1 > tol and o2 < -tol) or (o1 < -tol and o2 > tol)) and ((o3 > tol and o4 < -tol) or (o3 < -tol and o4 > tol)):
                raise StlGeometryError("extracted boundary self-intersects")


def _weld_vertices(raw: list[Point], tolerance: float) -> tuple[list[Point], list[tuple[int, int, int]]]:
    """Weld against first-seen representatives using exact Euclidean distance."""
    origin = tuple(min(p[axis] for p in raw) for axis in range(3))
    buckets: dict[tuple[int, int, int], list[int]] = defaultdict(list)
    points: list[Point] = []
    tri_ids: list[tuple[int, int, int]] = []
    for start in range(0, len(raw), 3):
        ids = []
        for p in raw[start:start+3]:
            relative = tuple((p[axis] - origin[axis]) / tolerance for axis in range(3))
            if not all(math.isfinite(value) for value in relative):
                raise StlGeometryError("coordinate range is too large for the requested weld tolerance")
            key = tuple(math.floor(value) for value in relative)
            candidates = []
            for dx in (-1, 0, 1):
                for dy in (-1, 0, 1):
                    for dz in (-1, 0, 1):
                        candidates.extend(buckets.get((key[0]+dx, key[1]+dy, key[2]+dz), ()))
            best = None
            best_distance = tolerance
            for point_id in candidates:
                q = points[point_id]
                distance = math.dist(p, q)
                if distance < best_distance or (
                    distance == best_distance and best is not None and point_id < best
                ):
                    best, best_distance = point_id, distance
            if best is None:
                best = len(points)
                points.append(p)
                buckets[key].append(best)
            ids.append(best)
        tri_ids.append(tuple(ids))
    return points, tri_ids

def extract_boundary(triangles: Iterable[Triangle], tolerance: float = 1e-9) -> list[Point]:
    if tolerance <= 0 or not math.isfinite(tolerance):
        raise StlGeometryError("tolerance must be positive and finite")
    raw = [p for tri in triangles for p in tri]
    if len(raw) < 3:
        raise StlGeometryError("STL contains no triangles")
    zs = [p[2] for p in raw]
    if max(zs)-min(zs) > tolerance:
        raise StlGeometryError("first-stage 2D STL must lie in a constant-Z XY plane")
    points, tri_ids = _weld_vertices(raw, tolerance)
    for ids in tri_ids:
        if len(set(ids)) != 3:
            raise StlGeometryError("degenerate triangle after vertex welding")
        a,b,c = (points[i] for i in ids)
        area2 = abs((b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0]))
        if area2 <= tolerance*tolerance:
            raise StlGeometryError("zero-area triangle in XY plane")
    incidence = Counter()
    for a,b,c in tri_ids:
        for u,v in ((a,b),(b,c),(c,a)):
            incidence[tuple(sorted((u,v)))] += 1
    if any(v > 2 for v in incidence.values()):
        raise StlGeometryError("non-manifold edge shared by more than two triangles")
    edges = [e for e,n in incidence.items() if n == 1]
    if len(edges) < 3:
        raise StlGeometryError("no closed outer boundary")
    adjacency = defaultdict(list)
    for a,b in edges:
        adjacency[a].append(b); adjacency[b].append(a)
    if any(len(v) != 2 for v in adjacency.values()):
        raise StlGeometryError("boundary is open, branched, or has touching loops")
    start = min(adjacency)
    order, previous, current = [start], None, start
    while True:
        candidates = [v for v in adjacency[current] if v != previous]
        nxt = min(candidates) if previous is None else candidates[0]
        if nxt == start: break
        if nxt in order: raise StlGeometryError("boundary traversal repeats a vertex")
        order.append(nxt); previous, current = current, nxt
        if len(order) > len(adjacency): raise StlGeometryError("boundary traversal did not close")
    if len(order) != len(adjacency):
        raise StlGeometryError("multiple boundary loops/components are not yet supported")
    loop = [points[i] for i in order]
    area2 = sum(loop[i][0]*loop[(i+1)%len(loop)][1]-loop[(i+1)%len(loop)][0]*loop[i][1] for i in range(len(loop)))
    if abs(area2) <= tolerance*tolerance: raise StlGeometryError("boundary encloses zero area")
    if area2 < 0: loop = [loop[0], *reversed(loop[1:])]
    _check_simple_loop(loop, tolerance)
    return loop

