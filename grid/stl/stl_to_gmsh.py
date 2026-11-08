"""Convert a planar triangulated STL region into Gmsh 2D geometry."""
from __future__ import annotations
import argparse, json, math, re, struct, subprocess
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable

Point = tuple[float, float, float]
Triangle = tuple[Point, Point, Point]

class Stl2DError(ValueError):
    """Unsupported or invalid planar STL input."""

def read_stl(path: Path) -> list[Triangle]:
    data = path.read_bytes()
    if len(data) >= 84:
        count = struct.unpack_from("<I", data, 80)[0]
        if 84 + count * 50 == len(data):
            result, offset = [], 84
            for _ in range(count):
                row = struct.unpack_from("<12fH", data, offset)
                result.append(tuple(tuple(float(v) for v in row[i:i+3]) for i in (3, 6, 9)))
                offset += 50
            if not all(math.isfinite(v) for tri in result for p in tri for v in p):
                raise Stl2DError("binary STL contains non-finite coordinates")
            return result
    try:
        content = data.decode("ascii")
    except UnicodeDecodeError as exc:
        raise Stl2DError("input is neither binary STL nor ASCII STL") from exc
    rows = re.findall(r"\bvertex\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)", content, re.I)
    if not rows or len(rows) % 3:
        raise Stl2DError("ASCII STL must contain complete groups of three vertices")
    points = [tuple(float(v) for v in row) for row in rows]
    if not all(math.isfinite(v) for p in points for v in p):
        raise Stl2DError("ASCII STL contains non-finite coordinates")
    return [tuple(points[i:i+3]) for i in range(0, len(points), 3)]

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
                raise Stl2DError("extracted boundary self-intersects")

def extract_boundary(triangles: Iterable[Triangle], tolerance: float = 1e-9) -> list[Point]:
    if tolerance <= 0 or not math.isfinite(tolerance):
        raise Stl2DError("tolerance must be positive and finite")
    raw = [p for tri in triangles for p in tri]
    if len(raw) < 3:
        raise Stl2DError("STL contains no triangles")
    zs = [p[2] for p in raw]
    if max(zs)-min(zs) > tolerance:
        raise Stl2DError("first-stage 2D STL must lie in a constant-Z XY plane")
    keys, points, tri_ids = {}, [], []
    for start in range(0, len(raw), 3):
        ids = []
        for p in raw[start:start+3]:
            key = tuple(round(v/tolerance) for v in p)
            if key not in keys:
                keys[key] = len(points)
                points.append(p)
            ids.append(keys[key])
        if len(set(ids)) != 3:
            raise Stl2DError("degenerate triangle after vertex welding")
        a,b,c = (points[i] for i in ids)
        area2 = abs((b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0]))
        if area2 <= tolerance*tolerance:
            raise Stl2DError("zero-area triangle in XY plane")
        tri_ids.append(tuple(ids))
    incidence = Counter()
    for a,b,c in tri_ids:
        for u,v in ((a,b),(b,c),(c,a)):
            incidence[tuple(sorted((u,v)))] += 1
    if any(v > 2 for v in incidence.values()):
        raise Stl2DError("non-manifold edge shared by more than two triangles")
    edges = [e for e,n in incidence.items() if n == 1]
    if len(edges) < 3:
        raise Stl2DError("no closed outer boundary")
    adjacency = defaultdict(list)
    for a,b in edges:
        adjacency[a].append(b); adjacency[b].append(a)
    if any(len(v) != 2 for v in adjacency.values()):
        raise Stl2DError("boundary is open, branched, or has touching loops")
    start = min(adjacency)
    order, previous, current = [start], None, start
    while True:
        candidates = [v for v in adjacency[current] if v != previous]
        nxt = min(candidates) if previous is None else candidates[0]
        if nxt == start: break
        if nxt in order: raise Stl2DError("boundary traversal repeats a vertex")
        order.append(nxt); previous, current = current, nxt
        if len(order) > len(adjacency): raise Stl2DError("boundary traversal did not close")
    if len(order) != len(adjacency):
        raise Stl2DError("multiple boundary loops/components are not yet supported")
    loop = [points[i] for i in order]
    area2 = sum(loop[i][0]*loop[(i+1)%len(loop)][1]-loop[(i+1)%len(loop)][0]*loop[i][1] for i in range(len(loop)))
    if abs(area2) <= tolerance*tolerance: raise Stl2DError("boundary encloses zero area")
    if area2 < 0: loop = [loop[0], *reversed(loop[1:])]
    _check_simple_loop(loop, tolerance)
    return loop

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
        if any(not name for name in names) or len(set(names)) != 4:
            raise Stl2DError("the four outer domain boundaries require four distinct non-empty names")
        body_name = str(config.get("body_boundary_name", "Wall")).strip()
        if not body_name or body_name in names:
            raise Stl2DError("body_boundary_name must be non-empty and distinct from outer boundary names")
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
    return config, assignments

def build_geo(loop: list[Point], assignments: dict[tuple[int,int],str], config: dict) -> str:
    size = float(config.get("mesh_size", 1.0))
    if not math.isfinite(size) or size <= 0:
        raise Stl2DError("mesh_size must be positive and finite")
    lines = ['// Generated from planar STL by OneFLOW.', 'SetFactory("Built-in");', f"lc = {size:.17g};"]
    if "_domain_bounds" in config:
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
            lines.append(f"Point({first_outer_point+offset}) = {{{x:.17g}, {y:.17g}, 0, lc}};")
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
        for side in ("left", "right", "bottom", "top"):
            name = re.sub(r"[^A-Za-z0-9_.-]+", "_", config["_outer_boundaries"][side])
            lines.append(f'Physical Curve("{name}") = {{{outer_line_ids[side]}}};')
    else:
        one_quad = bool(config.get("one_quad_cell", False))
        if one_quad and len(loop) != 4:
            raise Stl2DError("one_quad_cell requires exactly four boundary vertices")
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
    surface = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(config.get("surface_name","Fluid")))
    lines.append(f'Physical Surface("{surface}") = {{1}};')
    if config.get("one_quad_cell", False) and "_domain_bounds" not in config:
        lines.extend([f"Transfinite Curve {{{i}}} = 2;" for i in range(1,5)])
        lines.extend(["Transfinite Surface {1};", "Recombine Surface {1};"])
    lines.extend(["Mesh 2;", ""])
    return "\n".join(lines)

def convert(stl: Path, sidecar: Path, geo: Path, tolerance: float = 1e-9) -> dict:
    triangles = read_stl(stl)
    loop = extract_boundary(triangles, tolerance)
    config, assignments = _read_sidecar(sidecar, loop, tolerance)
    geo.parent.mkdir(parents=True, exist_ok=True)
    geo.write_text(build_geo(loop, assignments, config), encoding="utf-8", newline="\n")
    return {"triangles":len(triangles), "boundary_vertices":len(loop), "boundary_edges":len(assignments), "body_boundary_edges":len(loop) if "_domain_bounds" in config else 0, "geo":str(geo)}

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
            subprocess.run([args.gmsh,str(args.geo),"-2","-format","msh2","-o",str(args.msh)],check=True)
    except (OSError,json.JSONDecodeError,Stl2DError,subprocess.CalledProcessError) as exc:
        parser.exit(2,f"stl_to_gmsh: {exc}\n")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
