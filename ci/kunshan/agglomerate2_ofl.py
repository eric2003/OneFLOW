#!/usr/bin/env python3
"""Build the 2N rung of the M6 ladder by per-parent bisecting merges.

Feasibility is proven first (check_merge_dirs.py, mode 2N): a direction
assignment exists with zero violations over all 873,216 interior faces.  This
tool consumes that assignment and builds the coarse mesh:

  brick id          = 2 * parent + half   (half: 0 = low side of the cut)
  A(fine cell c)    = 2 * (c // 8) + half(c)      (parent*8 + octant in refine)

Faces, strictly from the refined mesh's own faces through A():
  * interior face with A(left) == A(right): interior to a brick -> dropped;
  * interior face with A(left) != A(right): kept, owner/neighbour from A();
  * boundary face: kept, BC identity inherited from that face.

No node-set matching, no BC fallback, no coordinate-based identity anywhere.
Ring order is rebuilt geometrically per emitted face (outward from the brick
that owns the face), and boundary faces are emitted first (the .ofl reader
convention).

Usage: agglomerate2_ofl.py <refined8.ofl> <out2N.ofl> <parent.ofl>

The bisect directions are solved on the PARENT (N) mesh, whose per-cell frames
are the constraint variables (this is exactly the network that
check_merge_dirs.py already verified PASS on M6).  The refined mesh is only
used for construction: fine cell -> parent is the refine provenance
(parent*8 + octant), and faces are consumed with real face IDs.

"""
import struct
import sys
from array import array

sys.path.insert(0, __file__.rsplit("/", 1)[0])
import refine_ofl as R
import check_merge_dirs as C


def build_assignment(z, mode="2N"):
    """Solve the merge-direction assignment on the PARENT mesh (N)."""
    cf = R.build_cell_faces(z)
    lat = [C.cell_lattice_positions(R.hex_corners(z, cf, c))
           for c in range(z.n_cells)]
    cons = []
    for f in range(z.n_bfaces, z.n_faces):
        l, r = z.left[f], z.right[f]
        got, err = C.analyse_face(R.face_ring(z, f), lat[l], lat[r])
        if got is None:
            R.fail("face %d unanalysable: %s" % (f, err))
        cons.append((f, l, r, C.allowed_matrix(*got, mode)))
    verdict, assign = C.solve_perm_network(z.n_cells, cons)
    if verdict != "PASS":
        R.fail("direction assignment failed: %s" % assign)
    return cons, assign


def main():
    mode = "2N"
    args = [a for a in sys.argv[1:] if not a.startswith("--mode")]
    if "--mode" in sys.argv:
        mode = sys.argv[sys.argv.index("--mode") + 1]
    if mode not in ("2N", "4N"):
        print("--mode must be 2N or 4N"); return 2
    src, dst = args[0], args[1]
    if len(args) < 3:
        print("usage: agglomerate2_ofl.py <refined8.ofl> <out.ofl> <parent.ofl> "
              "[--mode 2N|4N]")
        return 2
    z, pid, zone_type = R.parse_ofl(src)
    zp, _, _ = R.parse_ofl(sys.argv[3])
    if z.n_cells % 8:
        R.fail("cell count %d not a multiple of 8" % z.n_cells)
    P = z.n_cells // 8
    if zp.n_cells != P:
        R.fail("parent cell count %d != refined/8 %d" % (zp.n_cells, P))
    cons, d_parent = build_assignment(zp, mode)


    # --- A(c): fine cell -> brick -------------------------------------------
    # fine id = parent*8 + octant k (refine provenance).  The cut direction
    # d(parent) is a LATTICE axis; brick half = octant's coordinate on that
    # axis.  Child octant lattice coords come from CHILD_LATTICE.
    import agglomerate_ofl as A0
    CH = A0.CHILD_LATTICE
    A = array("i", bytes(4 * z.n_cells))     # per FINE cell -> merged cell
    if mode == "2N":
        for p in range(P):
            d = d_parent[p]
            for k in range(8):
                A[8 * p + k] = 2 * p + CH[k][d]
    else:
        # 4N: merge PAIRS along the lattice axis d -> 4 bricks per parent.
        # The pair index drops the d coordinate: the remaining two lattice
        # coords (in their original order) index the 2x2 plane.
        for p in range(P):
            d = d_parent[p]
            for k in range(8):
                a, b, c = CH[k]
                rest = [x for x in (a, b, c) if True]
                coords = [a, b, c]
                others = [coords[i] for i in range(3) if i != d]
                idx = 2 * others[0] + others[1]
                A[8 * p + k] = 4 * p + idx

    # --- faces ---------------------------------------------------------------
    # Coarse faces are built by MERGING fine faces into quads:
    #   * every parent face pf (refined faces 4*pf + i, i = 0..3) contributes
    #     two coarse quads, one per (brick-left, brick-right) group of its
    #     sub-faces;
    #   * every parent's brick interface (the four interior sheets between
    #     the two halves) contributes one coarse quad.
    # Group rule: corners are the nodes appearing EXACTLY ONCE across the
    # group's rings (shared midpoints appear twice and are not corners).
    # Expected counts: 2*nFaces_parent + nCells_parent faces,
    # 2*nB_parent boundary faces.
    from collections import defaultdict
    # refined layout (refine_ofl.py): ex faces first, 4 per PARENT face in
    # parent order [0, 4*nFaces_parent), then 12 interior faces per parent
    # cell.  The parent mesh supplies the counts; the refined mesh must agree.
    nFacesP = zp.n_faces
    if z.n_faces != 4 * nFacesP + 12 * zp.n_cells:
        R.fail("refined face count %d != 4*%d + 12*%d"
               % (z.n_faces, nFacesP, zp.n_cells))
    if z.n_bfaces != 4 * zp.n_bfaces:
        R.fail("refined boundary count %d != 4*%d" % (z.n_bfaces, zp.n_bfaces))

    groups = []           # (list_of_fine_faces, bc_type, bc_name)

    # 1) parent faces (both boundary and inter-parent), refined id 4*pf + i
    for pf in range(nFacesP):
        gm = defaultdict(list)
        for i in range(4):
            f = 4 * pf + i
            al, ar = A[z.left[f]], (A[z.right[f]] if f >= z.n_bfaces else -1)
            gm[(al, ar)].append(f)
        bnd = pf < zp.n_bfaces if hasattr(zp, "n_bfaces") else None
        for (al, ar), fs in gm.items():
            pass  # any size 1..4; the once-only corner rule handles merging
            bt = zp.bc_type[pf] if pf < zp.n_bfaces else None
            bn = zp.bc_name[pf] if pf < zp.n_bfaces else None
            groups.append((fs, bt, bn, al, ar))

    # 2) brick interfaces: interior sheets of the refined mesh.  The four
    # sheets of one interface alternate left/right orientation (refine stores
    # left = c*8+ia, right = c*8+ib), so the group key is the UNORDERED brick
    # pair; the emitted left/right orientation follows the first sheet.
    n_in_faces = z.n_faces - 4 * nFacesP
    gm = defaultdict(list)
    for f in range(4 * nFacesP, z.n_faces):
        al, ar = A[z.left[f]], A[z.right[f]]
        if al == ar:
            continue       # dropped: interior to a brick
        gm[(al, ar) if al < ar else (ar, al)].append(f)
    for key, fs in gm.items():
        # In 4N mode pair-pair boundaries consist of two sheets; in 2N mode
        # each brick interface consists of four sheets.  Any size 1..4 merges
        # into one coarse quad by the same boundary-walk rule.
        if not 1 <= len(fs) <= 4:
            R.fail("brick interface group has %d faces: %s" % (len(fs), fs[:8]))
        al = A[z.left[fs[0]]]
        ar = A[z.right[fs[0]]]
        groups.append((fs, None, None, al, ar))

    # --- merge each group into one face -------------------------------------
    # Ring order comes from the fine faces themselves: collect each fine ring's
    # directed edges, cancel interior edges (each appears once in both
    # directions), and walk the remaining boundary cycle.  Nodes that lie
    # collinearly between their two boundary neighbours (edge midpoints of the
    # parent) are dropped, leaving the true corner ring.  No angle sort, no
    # centroid, no assumption about the group size.
    def ring_from_group(fs):
        # Undirected edge multiplicity: interior edges of the group tiling
        # appear exactly twice (in either storage direction -- refine stores
        # the four interface sheets with the same chirality, so directed
        # cancellation is wrong); boundary edges appear exactly once.
        und = defaultdict(int)
        rings = [R.face_ring(z, f) for f in fs]
        for ring in rings:
            m = len(ring)
            for i in range(m):
                a, b = ring[i], ring[(i + 1) % m]
                und[(a, b) if a < b else (b, a)] += 1
        bnd = [e for e, c in und.items() if c == 1]
        if not bnd:
            R.fail("group %s has no boundary edges" % (fs,))
        for e, c in und.items():
            if c > 2:
                R.fail("group %s edge %s multiplicity %d" % (fs, e, c))
        # adjacency walk: each boundary node has exactly two boundary edges
        adj = defaultdict(list)
        for (a, b) in bnd:
            adj[a].append(b)
            adj[b].append(a)
        for n, ns in adj.items():
            if len(ns) != 2:
                R.fail("group %s boundary node %s has %d edges"
                       % (fs, n, len(ns)))
        start = min(adj)
        walk = [start]
        prev = None
        cur = start
        while True:
            nxt = adj[cur][0] if adj[cur][0] != prev else adj[cur][1]
            if nxt == start:
                break
            walk.append(nxt)
            prev, cur = cur, nxt
        if len(walk) != len(adj):
            R.fail("group %s walk %d != nodes %d" % (fs, len(walk), len(adj)))
        # Corners by multiplicity: nodes of the boundary walk that appear in
        # exactly ONE fine face of the group are the merged face's corners
        # (edge midpoints and interior centre nodes appear twice or more).
        # The walk supplies the cyclic order.
        cnt = defaultdict(int)
        for ring in rings:
            for n in ring:
                cnt[n] += 1
        corners = [n for n in walk if cnt[n] == 1]
        if len(corners) != 4:
            R.fail("group %s has %d corners (walk %s, counts %s)"
                   % (fs, len(corners), walk, dict(cnt)))
        # order the corners by their position in the walk
        pos = {n: i for i, n in enumerate(walk)}
        corners.sort(key=lambda n: pos[n])
        return corners

    # emit: boundary groups first
    face_rings = []
    out_left = array("i")
    out_right = array("i")
    out_bc_type = array("i")
    out_bc_name = array("i")
    boundary_groups = [g for g in groups if g[1] is not None]
    interior_groups = [g for g in groups if g[1] is None]
    # expected boundary count from the assignment: a parent boundary face
    # contributes 1 coarse face when the bisect leaves it uncut (d == face
    # normal) and 2 when cut (d != normal).  Build the expectation from the
    # same A() used for emission.
    expected_b = 0
    for pf in range(zp.n_bfaces):
        f = 4 * pf                      # any sub-face identifies the pair
        seen = set()
        for i in range(4):
            seen.add(A[z.left[4 * pf + i]])
        expected_b += len(seen)
    if len(boundary_groups) != expected_b:
        R.fail("boundary groups %d != assignment-derived %d"
               % (len(boundary_groups), expected_b))
    for fs, bt, bn, al, ar in boundary_groups:
        face_rings.append(ring_from_group(fs))
        out_left.append(al)
        out_right.append(-1)
        out_bc_type.append(bt)
        out_bc_name.append(bn)
    n_b_out = len(out_bc_type)
    for fs, bt, bn, al, ar in interior_groups:
        face_rings.append(ring_from_group(fs))
        out_left.append(al)
        out_right.append(ar)

    # --- orient every emitted face outward from its left brick ---------------
    n_bricks = (2 if mode == "2N" else 4) * P
    cf = R.build_cell_faces(z)
    bcx = array("d", bytes(8 * n_bricks))
    bcy = array("d", bytes(8 * n_bricks))
    bcz = array("d", bytes(8 * n_bricks))
    bcnt = array("i", bytes(4 * n_bricks))
    for c in range(z.n_cells):
        nodes = set()
        for f in cf[c]:
            nodes.update(R.face_ring(z, f))
        b = A[c]
        for n in nodes:
            bcx[b] += z.xn[n]; bcy[b] += z.yn[n]; bcz[b] += z.zn[n]
        bcnt[b] += len(nodes)
    for b in range(n_bricks):
        if bcnt[b]:
            bcx[b] /= bcnt[b]; bcy[b] /= bcnt[b]; bcz[b] /= bcnt[b]

    # orientation: every emitted face's ring must point from its left brick to
    # its right brick (or outward at a boundary)
    for i in range(len(face_rings)):
        ring = face_rings[i]
        m = len(ring)
        fx = sum(z.xn[n] for n in ring) / m
        fy = sum(z.yn[n] for n in ring) / m
        fz = sum(z.zn[n] for n in ring) / m
        # normal by Newell's method (works for any planar polygon ring)
        nx_ = ny_ = nz_ = 0.0
        for k in range(m):
            a, b = ring[k], ring[(k + 1) % m]
            nx_ += (z.yn[a] - z.yn[b]) * (z.zn[a] + z.zn[b])
            ny_ += (z.zn[a] - z.zn[b]) * (z.xn[a] + z.xn[b])
            nz_ += (z.xn[a] - z.xn[b]) * (z.yn[a] + z.yn[b])
        l = out_left[i]
        r = out_right[i]
        if r >= 0:
            ox = bcx[r] - bcx[l]
            oy = bcy[r] - bcy[l]
            oz = bcz[r] - bcz[l]
        else:
            ox = fx - bcx[l]
            oy = fy - bcy[l]
            oz = fz - bcz[l]
        if nx_ * ox + ny_ * oy + nz_ * oz < 0.0:
            face_rings[i] = list(reversed(ring))

    new = R.Zone()
    new.n_nodes = z.n_nodes
    new.n_faces = len(face_rings)
    new.n_cells = n_bricks
    new.vol_bc_type = z.vol_bc_type
    new.xn, new.yn, new.zn = z.xn, z.yn, z.zn
    new.num_face_node = array("i", [len(r) for r in face_rings])
    new.fnmem = array("i")
    for r in face_rings:
        new.fnmem.extend(r)
    new.left = out_left
    new.right = out_right
    new.n_bfaces = n_b_out
    new.bc_type = out_bc_type
    new.bc_name = out_bc_name

    print("2N mesh: cells=%d faces=%d boundary=%d (groups=%d)"
          % (new.n_cells, new.n_faces, new.n_bfaces, len(groups)))
    R.write_ofl(dst, new, pid, zone_type)
    return 0


if __name__ == "__main__":
    sys.exit(main())
