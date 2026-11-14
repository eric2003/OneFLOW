#!/usr/bin/env python3
"""Agglomerate a uniformly refined .ofl mesh back to a coarser ladder rung.

The uniform 2x2x2 refinement numbers children as parent*8 + k, where k is the
parent frame corner index:

    0=(0,0,0) 1=(1,0,0) 2=(1,1,0) 3=(0,1,0)
    4=(0,0,1) 5=(1,0,1) 6=(1,1,1) 7=(0,1,1)        (lattice coords a,b,c)

so every child carries its own lattice position and no direction labelling is
needed.  Merging groups of children back into bricks therefore inherits the
input's topology exactly - including the M6 O-grid pole folds that break a
forward one-direction cut.

Group size 2 (pattern full x full x half, i.e. one direction merged) halves the
cell count; groups of 4 (full x half x half) quarter it.  With N = 294,912:

    merge 8N with group 2 -> 4N   (1,179,648 cells)
    merge 8N with group 4 -> 2N   (  589,824 cells)
    merge 64N with group 2 -> 32N (9,437,184 cells)
    merge 64N with group 4 -> 16N (4,718,592 cells)

Lattice node recovery: the node at lattice (a,b,c) is the only node shared by
every child whose (i,j,k) range contains it, so it is the intersection of
those children's corner sets.  Faces are then deduplicated by node set and
oriented outward from the owning brick's centroid; boundary sub-faces inherit
bcType/bcNameId from the parent boundary face that contains them.

Usage: python3 agglomerate_ofl.py <refined.ofl> <output.ofl> --group 2|4
"""
import struct
import sys
from array import array

import refine_ofl as R

CHILD_LATTICE = ((0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
                 (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1))


def main():
    if len(sys.argv) < 4:
        print(__doc__)
        return 2
    src, dst = sys.argv[1], sys.argv[2]
    group = 0
    for i in range(3, len(sys.argv)):
        if sys.argv[i] == "--group":
            group = int(sys.argv[i + 1])
    if group not in (2, 4):
        print("--group must be 2 or 4")
        return 2

    report = []
    z, pid, zone_type = R.parse_ofl(src)
    if z.n_cells % 8 != 0:
        R.fail("input cell count %d is not a multiple of 8" % z.n_cells)
    parents = z.n_cells // 8
    _, vol_in = R.cell_volumes(z, report)

    cell_faces = R.build_cell_faces(z)
    corners = array("i", bytes(4 * 8 * z.n_cells))
    for c in range(z.n_cells):
        corners[8 * c:8 * c + 8] = array("i", R.hex_corners(z, cell_faces, c))

    # brick pattern in child units: group 2 -> (2,2,1); group 4 -> (2,1,1)
    dims = (2, 2, 1) if group == 2 else (2, 1, 1)

    def lattice_map(pset):
        """Map lattice (a,b,c) in {0,1,2}^3 to the node id inside one parent.

        Every node of the parent's 3x3x3 lattice is shared by a known set of
        children (a corner by 1, an edge midpoint by 2, a face centre by 4, the
        cell centre by 8), and its lattice coordinate follows from whether that
        set straddles the mid-plane in each axis: all-zero -> 0, all-one -> 2,
        mixed -> 1.  Purely combinatorial, no geometry, so collapsed O-grid
        layers cannot confuse it.
        """
        owners = {}
        for k in range(8):
            for n in pset[k]:
                owners.setdefault(n, []).append(k)
        out = {}
        for n, ks in owners.items():
            coord = []
            for axis in range(3):
                vals = {CHILD_LATTICE[k][axis] for k in ks}
                if vals == {0}:
                    coord.append(0)
                elif vals == {1}:
                    coord.append(2)
                else:
                    coord.append(1)
            out[tuple(coord)] = n
        if len(out) != 27:
            R.fail("parent lattice produced %d positions, expected 27" % len(out))
        return out

    # Per-parent boundary records: at a collapsed cell a merged sub-face need
    # not be a node subset of one parent face, so the parent's own boundary
    # faces are the reliable source.
    pbc = {}
    for f in range(z.n_bfaces):
        pbc.setdefault(z.left[f] // 8, []).append(
            (frozenset(R.face_ring(z, f)), z.bc_type[f], z.bc_name[f]))

    # parent boundary faces for bc inheritance
    bnode = {}
    for f in range(z.n_bfaces):
        ring = R.face_ring(z, f)
        key = frozenset(ring)
        for n in ring:
            bnode.setdefault(n, []).append((key, z.bc_type[f], z.bc_name[f]))

    def cache_ensure(cache, pset, lookup):
        if lookup not in cache:
            if "__lat__" not in cache:
                cache["__lat__"] = lattice_map(pset)
            cache[lookup] = cache["__lat__"][lookup]
        return cache[lookup]

    brick_corners = []
    brick_centroid = array("d")
    face_cells = {}
    for p in range(parents):
        pset = [set(corners[8 * (8 * p + k):8 * (8 * p + k) + 8]) for k in range(8)]
        cache = {}
        for a0 in range(0, 2, dims[0]):
            for b0 in range(0, 2, dims[1]):
                for c0 in range(0, 2, dims[2]):
                    bid = len(brick_corners)
                    cs = []
                    for a in (a0, min(a0 + dims[0], 2)):
                        for b in (b0, min(b0 + dims[1], 2)):
                            for c in (c0, min(c0 + dims[2], 2)):
                                cs.append(cache_ensure(cache, pset, (a, b, c)))
                    brick_corners.append(cs)
                    sx = sum(z.xn[n] for n in cs) / 8.0
                    sy = sum(z.yn[n] for n in cs) / 8.0
                    sz = sum(z.zn[n] for n in cs) / 8.0
                    brick_centroid.extend((sx, sy, sz))
                    # six faces of a unit lattice box: (a,b,c) axis and offset
                    for axis in range(3):
                        for side in (0, 1):
                            quad = []
                            for (u, v) in ((0, 0), (1, 0), (1, 1), (0, 1)):
                                coord = [None, None, None]
                                coord[axis] = side
                                others = [x for x in range(3) if x != axis]
                                coord[others[0]] = u
                                coord[others[1]] = v
                                lookup = (a0 + coord[0] * dims[0],
                                          b0 + coord[1] * dims[1],
                                          c0 + coord[2] * dims[2])
                                quad.append(cache_ensure(cache, pset, lookup))
                            key = frozenset(quad)
                            face_cells.setdefault(key, []).append((bid, quad))
    # emit
    out_frings = array("i")
    out_left = array("i")
    out_right = array("i")
    out_bc_type = array("i")
    out_bc_name = array("i")
    n_boundary = 0
    pairs = []
    for key, owners in face_cells.items():
        # Collapsed cells make several lattice faces share one node set, so a
        # node-set bucket can hold more than two entries and can hold two faces
        # of the SAME brick.  Only faces from different bricks form an
        # interface; whatever is left is a boundary face of its own brick.
        rest = list(owners)
        used = [False] * len(rest)
        for i in range(len(rest)):
            if used[i]:
                continue
            paired = False
            for j in range(i + 1, len(rest)):
                if used[j] or rest[j][0] == rest[i][0]:
                    continue
                pairs.append((rest[i][0], rest[i][1], rest[j][0]))
                used[i] = used[j] = True
                paired = True
                break
            if not paired:
                # No partner from a different brick: a boundary face, or one of
                # several coincident faces of a collapsed brick (which cancel in
                # the volume exactly like the input's fold duplicates).
                pairs.append((rest[i][0], rest[i][1], None))
                used[i] = True
    # Second pass: pair leftovers through a node bucket (linear, not O(L^2)).
    # A genuine interface face shares all four nodes; three still distinguishes
    # it from a neighbour that merely touches an edge.
    node_bucket = {}
    for idx, pr in enumerate(pairs):
        if pr[2] is None:
            for n in pr[1]:
                node_bucket.setdefault(n, []).append(idx)
    taken = set()
    n_repair = 0
    for idx, pr in enumerate(pairs):
        if pr[2] is not None or idx in taken:
            continue
        qi = set(pr[1])
        partner = None
        for n in pr[1]:
            for j in node_bucket.get(n, ()):
                if j <= idx or j in taken or pairs[j][0] == pr[0]:
                    continue
                if len(qi & set(pairs[j][1])) >= 3:
                    partner = j
                    break
            if partner is not None:
                break
        if partner is not None:
            pairs[idx] = (pr[0], pr[1], pairs[partner][0])
            taken.add(idx)
            taken.add(partner)
            n_repair += 1
    report.append("second-pass interface repair: %d face(s)" % n_repair)

    for bid, quad, other in pairs:
        quad_set = frozenset(quad)
        # orient: normal from `bid` toward the neighbour (or outward)
        fx = sum(z.xn[n] for n in quad) / 4.0
        fy = sum(z.yn[n] for n in quad) / 4.0
        fz = sum(z.zn[n] for n in quad) / 4.0
        d1x = z.xn[quad[2]] - z.xn[quad[0]]
        d1y = z.yn[quad[2]] - z.yn[quad[0]]
        d1z = z.zn[quad[2]] - z.zn[quad[0]]
        d2x = z.xn[quad[3]] - z.xn[quad[1]]
        d2y = z.yn[quad[3]] - z.yn[quad[1]]
        d2z = z.zn[quad[3]] - z.zn[quad[1]]
        nx_ = d1y * d2z - d1z * d2y
        ny_ = d1z * d2x - d1x * d2z
        nz_ = d1x * d2y - d1y * d2x
        if other is not None:
            ox = brick_centroid[3 * other] - brick_centroid[3 * bid]
            oy = brick_centroid[3 * other + 1] - brick_centroid[3 * bid + 1]
            oz = brick_centroid[3 * other + 2] - brick_centroid[3 * bid + 2]
        else:
            ox = fx - brick_centroid[3 * bid]
            oy = fy - brick_centroid[3 * bid + 1]
            oz = fz - brick_centroid[3 * bid + 2]
        ring = quad if nx_ * ox + ny_ * oy + nz_ * oz >= 0.0 else \
            (quad[0], quad[3], quad[2], quad[1])
        out_frings.extend(ring)
        out_left.append(bid)
        if other is not None:
            out_right.append(other)
        else:
            n_boundary += 1
            out_right.append(-1)
            bc = None
            cands = []
            for n in quad:
                cands.extend(bnode.get(n, ()))
            for (fset, bt, bn) in cands:
                if quad_set <= fset:
                    bc = (bt, bn)
                    break
            if bc is None:
                # Degenerate parents collapse lattice nodes, so a sub-face may
                # not be a strict subset of one parent face; take the parent
                # boundary face sharing the most nodes instead.
                best = 0
                for (fset, bt, bn) in cands:
                    overlap = len(quad_set & fset)
                    if overlap > best:
                        best = overlap
                        bc = (bt, bn)
                if bc is None or best < 2:
                    parent = bid // group
                    for (fset, bt, bn) in pbc.get(parent, ()):
                        overlap = len(quad_set & fset)
                        if overlap > best:
                            best = overlap
                            bc = (bt, bn)
                if bc is None:
                    R.fail("merged boundary face has no parent bc record")
                report.append("bc fallback: subset failed, overlap=%d" % best)
            out_bc_type.append(bc[0])
            out_bc_name.append(bc[1])

    new = R.Zone()
    new.n_nodes = z.n_nodes
    new.n_faces = len(out_frings) // 4
    new.n_cells = len(brick_corners)
    new.vol_bc_type = z.vol_bc_type
    new.xn, new.yn, new.zn = z.xn, z.yn, z.zn
    new.num_face_node = array("i", b"\x04\x00\x00\x00" * new.n_faces)
    new.fnmem = out_frings
    new.left = out_left
    new.right = out_right
    new.n_bfaces = n_boundary
    new.bc_type = out_bc_type
    new.bc_name = out_bc_name

    report.append("agglomerate: parents=%d group=%d -> cells=%d boundary=%d"
                  % (parents, group, new.n_cells, n_boundary))
    if new.n_cells != z.n_cells // group:
        R.fail("cell count %d != %d" % (new.n_cells, z.n_cells // group))
    R.write_ofl(dst, new, pid, zone_type)
    R.validate_output(dst, vol_in, report)
    for line in report:
        print(line)
    return 0


if __name__ == "__main__":
    sys.exit(main())
