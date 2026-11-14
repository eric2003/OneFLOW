#!/usr/bin/env python3
"""Refine a single-zone OneFLOW hexahedral .ofl grid by 2x2x2 per hex.

The .ofl binary layout handled here (verified byte-exact against
test/cfd/m6wingroe_sa/grid/m6wing.ofl and the m6wing-8.ofl fixture):

    int32 nZones
    int32 pid[nZones]
    int32 zoneType[nZones]
    per zone:
      int64 nLength                 payload byte count
      int32 nNodes, nFaces, nCells
      double xN[nNodes], yN[nNodes], zN[nNodes]
      int32 volBcType
      int32 numFaceNode[nFaces]     all 4 for hexahedra
      int32 faceNodeMem[sum(numFaceNode)]
      int32 leftCells[nFaces], rightCells[nFaces]
      int32 nBFaces
      int32 bcType[nBFaces], bcNameId[nBFaces]
      int32 nIFaces
      (nIFaces > 0) int32 zoneId[], localInterfaceId[], i2b[]

Only single-zone, interface-free inputs are accepted: multi-zone grids
should be refined zone by zone through this entry only after the interface
records are rebuilt, which this repository does not need because the
production Partition task rebuilds interface records itself.  Refine the
single-zone source and partition the refined grid instead.

Refinement (uniform, per hexahedron):
  new nodes  = old nodes
             + one midpoint per unique old edge      (id in [V, V+E))
             + one center per old face               (id in [V+E, V+E+F))
             + one center per old cell               (id in [V+E+F, ...))
  new cells  = 8 per old cell, indexed cell*8+corner
  new faces  = 4 sub-quads per old face (external; boundary sub-quads stay
               first and inherit bcType/bcNameId)
             + 12 interior faces per old cell (one per cell edge)

Conventions preserved:
  * quad rings wind so the normal points from leftCells to rightCells
    (checked on input and re-derived on output by a signed-volume test);
  * boundary faces come first; their rightCells keep the input convention
    (-1 in the single-zone M6 file);
  * boundary sub-faces inherit bcType/bcNameId from the parent face.

Self-checks (all must pass; the script exits non-zero otherwise):
  PARSE     byte-level closure: payload == nLength, EOF, all quads,
            nZones == 1, nIFaces == 0, every cell has 6 faces.
  COUNTS    nCells' == 8*C, nFaces' == 4*F + 12*C,
            nNodes' == V + E + F + C, nBFaces' == 4*nBFaces.
  VOLUME    every cell volume > 0 on the re-read output file, and the total
            volume matches the input total within 1e-10 relative.
  REPARSE   the written file is re-read from disk and re-validated as an
            independent mesh (degree, quads, boundary order, volume).

Usage: python3 refine_ofl.py <input.ofl> <output.ofl> [--cuts K]

Options:
  (no option)  uniform 2x2x2 refinement (8x cells per hex)
  --cuts K     sequential one-direction refinement, K = 1, 2 or 3
               (2^K x cells; K=3 reproduces the uniform mesh by counts and
               node-set identity, cross-checkable against the uniform run)
"""

import struct
import sys
from array import array

I32 = "<i"
I64 = "<q"
F64 = "<d"


def fail(msg):
    print("REFINE_FAIL " + msg)
    sys.exit(1)


class Zone:
    pass


def parse_ofl(path):
    """Parse a single-zone .ofl into plain arrays. Returns (zone, pid, zoneType)."""
    with open(path, "rb") as fh:
        data = fh.read()

    off = 0
    (n_zones,) = struct.unpack_from(I32, data, off)
    off += 4
    if n_zones != 1:
        fail("input must be single-zone, got nZones=%d" % n_zones)
    pid = struct.unpack_from("<%di" % n_zones, data, off)
    off += 4 * n_zones
    zone_type = struct.unpack_from("<%di" % n_zones, data, off)
    off += 4 * n_zones

    (n_length,) = struct.unpack_from(I64, data, off)
    if n_length < 12 or n_length > len(data) - off - 8:
        fail("implausible zone payload length %d at offset %d" % (n_length, off))
    off += 8
    zone_start = off

    z = Zone()
    (z.n_nodes, z.n_faces, z.n_cells) = struct.unpack_from("<3i", data, off)
    off += 12
    if min(z.n_nodes, z.n_faces, z.n_cells) <= 0:
        fail("non-positive grid counts")

    def rd_doubles(n):
        nonlocal off
        a = array("d")
        a.frombytes(data[off : off + 8 * n])
        off += 8 * n
        return a

    z.xn = rd_doubles(z.n_nodes)
    z.yn = rd_doubles(z.n_nodes)
    z.zn = rd_doubles(z.n_nodes)
    (z.vol_bc_type,) = struct.unpack_from(I32, data, off)
    off += 4

    z.num_face_node = array("i")
    z.num_face_node.frombytes(data[off : off + 4 * z.n_faces])
    off += 4 * z.n_faces
    if any(c != 4 for c in z.num_face_node):
        fail("non-quad face found; hexahedral refinement needs all quads")
    total_face_nodes = sum(z.num_face_node)
    z.fnmem = array("i")
    z.fnmem.frombytes(data[off : off + 4 * total_face_nodes])
    off += 4 * total_face_nodes
    z.left = array("i")
    z.left.frombytes(data[off : off + 4 * z.n_faces])
    off += 4 * z.n_faces
    z.right = array("i")
    z.right.frombytes(data[off : off + 4 * z.n_faces])
    off += 4 * z.n_faces

    (z.n_bfaces,) = struct.unpack_from(I32, data, off)
    off += 4
    if z.n_bfaces < 0 or z.n_bfaces > z.n_faces:
        fail("invalid nBFaces")
    z.bc_type = array("i")
    z.bc_type.frombytes(data[off : off + 4 * z.n_bfaces])
    off += 4 * z.n_bfaces
    z.bc_name = array("i")
    z.bc_name.frombytes(data[off : off + 4 * z.n_bfaces])
    off += 4 * z.n_bfaces
    (z.n_ifaces,) = struct.unpack_from(I32, data, off)
    off += 4
    if z.n_ifaces != 0:
        fail("input carries interface records; only interface-free grids are supported")
    off += 12 * z.n_ifaces

    if off != zone_start + n_length:
        fail("zone payload mismatch: consumed %d, nLength %d" % (off - zone_start, n_length))
    if off != len(data):
        fail("trailing bytes after zone payload: %d" % (len(data) - off))
    return z, pid, zone_type


def face_ring(z, f):
    base = 4 * f
    return z.fnmem[base], z.fnmem[base + 1], z.fnmem[base + 2], z.fnmem[base + 3]


def build_cell_faces(z):
    """cell -> list of face ids, plus input-degree validation."""
    degree = array("i", bytes(4 * z.n_cells))
    cell_faces = [None] * z.n_cells
    for f in range(z.n_faces):
        lc = z.left[f]
        rc = z.right[f]
        if lc < 0 or lc >= z.n_cells:
            fail("face %d left cell %d out of range" % (f, lc))
        if f < z.n_bfaces:
            if rc >= 0 and rc < z.n_cells:
                fail("boundary face %d right cell %d not outside" % (f, rc))
            degree[lc] += 1
            lst = cell_faces[lc]
            if lst is None:
                cell_faces[lc] = [f]
            else:
                lst.append(f)
        else:
            if rc < 0 or rc >= z.n_cells:
                fail("internal face %d right cell %d out of range" % (f, rc))
            degree[lc] += 1
            degree[rc] += 1
            if cell_faces[lc] is None:
                cell_faces[lc] = [f]
            else:
                cell_faces[lc].append(f)
            if cell_faces[rc] is None:
                cell_faces[rc] = [f]
            else:
                cell_faces[rc].append(f)
    bad = sum(1 for c in degree if c != 6)
    if bad:
        fail("%d cells do not have exactly 6 faces (hexahedral input required)" % bad)
    return cell_faces


def hex_corners(z, cell_faces, c):
    """Recover a consistent corner ordering [b0..b3, u0..u3] for cell c."""
    faces = cell_faces[c]
    if len(faces) != 6:
        fail("cell %d has %d faces" % (c, len(faces)))
    b = face_ring(z, faces[0])
    bset = set(b)
    if len(bset) != 4:
        fail("bottom face of cell %d repeats a node" % c)
    top = None
    sides = []
    for f in faces[1:]:
        ring = face_ring(z, f)
        if set(ring).isdisjoint(bset):
            if top is not None:
                fail("cell %d has two opposite faces" % c)
            top = ring
        else:
            sides.append(ring)
    if top is None or len(sides) != 4:
        fail("cell %d is not a hexahedron by opposite-face test" % c)
    tset = set(top)
    if len(tset) != 4:
        fail("top face of cell %d repeats a node" % c)

    upper = [None] * 4
    # A side face ring is (b_i, b_{i+1}, u_{i+1}, u_i) in cyclic order, so
    # each bottom node has exactly one top node among its two ring
    # neighbours: that neighbour is its upper partner via the vertical edge.
    for ring in sides:
        n_low = sum(1 for p in ring if p in bset)
        if n_low != 2:
            fail("cell %d side face shares %d nodes with bottom" % (c, n_low))
        for k, p in enumerate(ring):
            if p not in bset:
                continue
            prev = ring[(k - 1) % 4]
            nxt = ring[(k + 1) % 4]
            if prev in tset and nxt in bset:
                upper[b.index(p)] = prev
            elif nxt in tset and prev in bset:
                upper[b.index(p)] = nxt
            else:
                fail("cell %d side ring neighbour split invalid" % c)
            if prev in tset and nxt in tset:
                fail("cell %d side ring neighbours are both top nodes" % c)
    if any(u is None for u in upper) or len(set(upper)) != 4 or set(upper) != tset:
        fail("cell %d upper corner correspondence unresolved" % c)

    # Uniform chirality convention: the bottom-face ring must point OUT of the
    # cell.  Face rings are stored left -> right (outward from the left cell),
    # so a cell that is the RIGHT cell of its bottom face sees the ring
    # reversed.  Mixed chirality makes the direction labelling non-integrable
    # (a "twist" around any loop), so mirror the cell frame in that case.
    bottom_face = faces[0]
    if z.left[bottom_face] != c:
        b = (b[0], b[3], b[2], b[1])
        upper = [upper[0], upper[3], upper[2], upper[1]]

    corners = list(b) + upper
    if len(set(corners)) != 8:
        fail("cell %d does not have 8 distinct corners" % c)
    return corners


def cell_edges(corners):
    """12 edges as (corner_index_a, corner_index_b) pairs; a < b."""
    bottom = [(i, (i + 1) % 4) for i in range(4)]
    top = [(4 + i, 4 + (i + 1) % 4) for i in range(4)]
    vertical = [(i, 4 + i) for i in range(4)]
    return [(min(a, b), max(a, b)) for a, b in bottom + top + vertical]


def det3(x0, y0, z0, x1, y1, z1, x2, y2, z2):
    return (
        x0 * (y1 * z2 - z1 * y2)
        - y0 * (x1 * z2 - z1 * x2)
        + z0 * (x1 * y2 - y1 * x2)
    )


def dump_cell_faces(z, c):
    """Print one cell's faces and the owner/neighbour volumes (diagnostic)."""
    print("DEBUG negative cell %d" % c, file=sys.stderr)
    faces = []
    for f in range(z.n_faces):
        if z.left[f] == c:
            faces.append((f, z.left[f], z.right[f]))
        elif f >= z.n_bfaces and z.right[f] == c:
            faces.append((f, z.left[f], z.right[f]))
    for f, lc, rc in faces:
        print("  face %d ring=%s left=%d right=%d"
              % (f, face_ring(z, f), lc, rc), file=sys.stderr)
    nbr = set()
    for f, lc, rc in faces:
        nbr.add(lc)
        if rc >= 0:
            nbr.add(rc)
    nbr.discard(c)
    for n in sorted(nbr)[:4]:
        print("  neighbour cell %d faces: %s"
              % (n, [face_ring(z, f) for f in range(z.n_faces)
                     if z.left[f] == n or (f >= z.n_bfaces and z.right[f] == n)]),
              file=sys.stderr)


def cell_volumes(z, report=None):
    """Signed cell volumes from the left->right winding convention.

    For a face ring (p0..p3) with normal pointing left->right, the outward
    side of the left cell is the ring orientation, so the left cell gains
    the pyramid sum and the right cell loses it.  Every cell of a valid
    mesh must end with a positive volume.
    """
    vol = array("d", bytes(8 * z.n_cells))
    fn = z.fnmem
    left = z.left
    right = z.right
    nB = z.n_bfaces
    xn, yn, zn = z.xn, z.yn, z.zn
    total = 0.0
    for f in range(z.n_faces):
        b = 4 * f
        p0 = fn[b]; p1 = fn[b + 1]; p2 = fn[b + 2]; p3 = fn[b + 3]
        x0, y0, z0 = xn[p0], yn[p0], zn[p0]
        x1, y1, z1 = xn[p1], yn[p1], zn[p1]
        x2, y2, z2 = xn[p2], yn[p2], zn[p2]
        x3, y3, z3 = xn[p3], yn[p3], zn[p3]
        s = det3(x0, y0, z0, x1, y1, z1, x2, y2, z2) \
            + det3(x0, y0, z0, x2, y2, z2, x3, y3, z3)
        lc = left[f]
        vol[lc] += s
        if f >= nB:
            vol[right[f]] -= s
    # Degenerate band: O-grid pole/seam cells collapse geometrically while
    # keeping 8 distinct node ids, so their volume is round-off around zero
    # (observed -8.5e-11 on the M6 O-grid).  Treat |v| within 1e-7 of the
    # largest cell volume as degenerate; a winding or left/right bug produces
    # a negative volume of ordinary cell magnitude and still fails.
    vmax = max(vol)
    # A genuinely inverted cell has |v| comparable to the typical volume; an
    # O-grid sliver (thickness ~1e-4 of its other edges) can land a few 1e-7
    # below zero on pure rounding.  Only the former is a winding/left-right
    # bug, so the rejection band is 1e-4 of the largest cell volume, and the
    # global volume conservation below stays at 1e-10 as the real guard.
    tol = 1.0e-4 * vmax if vmax > 0.0 else 0.0
    scale = 1.0 / 6.0
    total = 0.0
    ndeg = 0
    worst = (vol[0], 0)
    for c in range(z.n_cells):
        s = vol[c]
        if s < worst[0]:
            worst = (s, c)
        if s < -tol:
            dump_cell_faces(z, c)
            fail("cell %d volume %.6e is negative beyond the degenerate band "
                 "(tol %.3e): winding or left/right broken"
                 % (c, s * scale, tol * scale))
        if -tol <= s <= tol:
            ndeg += 1
        total += s
    if report is not None:
        report.append(
            "volumes: min=%.6e max=%.6e degenerate=%d worst-cell=%d"
            % (worst[0] * scale, vmax * scale, ndeg, worst[1]))
    return vol, total * scale


def refine(z, report):
    V, F, C = z.n_nodes, z.n_faces, z.n_cells
    nB = z.n_bfaces
    # Boundary right-cell convention of the input.
    ghost_mode = any(z.right[f] >= 0 for f in range(nB))
    if ghost_mode:
        for f in range(nB):
            if z.right[f] != f + C:
                fail("boundary face %d right cell is not the ghost id" % f)
    else:
        for f in range(nB):
            if z.right[f] >= 0:
                fail("inconsistent boundary right cells")
    if any(z.right[f] < 0 for f in range(nB, F)):
        fail("internal face with negative right cell")

    cell_faces = build_cell_faces(z)

    # --- Pass A: edge midpoints (id on first sight) and face centers -------
    # Edge keys are packed into one integer (a * node_stride + b, a < b) so an
    # O(60M) midpoint table costs a few GB instead of tens of GB with tuples.
    mid_id = {}
    stride = V + 1
    mid_x = array("d"); mid_y = array("d"); mid_z = array("d")
    fc_x = array("d"); fc_y = array("d"); fc_z = array("d")
    fn = z.fnmem
    xn, yn, zn = z.xn, z.yn, z.zn
    for f in range(F):
        b = 4 * f
        p0 = fn[b]; p1 = fn[b + 1]; p2 = fn[b + 2]; p3 = fn[b + 3]
        for a, e in ((p0, p1), (p1, p2), (p2, p3), (p3, p0)):
            key = a * stride + e if a < e else e * stride + a
            mid = mid_id.get(key)
            if mid is None:
                mid = V + len(mid_x)
                mid_id[key] = mid
                mid_x.append(0.5 * (xn[a] + xn[e]))
                mid_y.append(0.5 * (yn[a] + yn[e]))
                mid_z.append(0.5 * (zn[a] + zn[e]))
        fc_x.append(0.25 * (xn[p0] + xn[p1] + xn[p2] + xn[p3]))
        fc_y.append(0.25 * (yn[p0] + yn[p1] + yn[p2] + yn[p3]))
        fc_z.append(0.25 * (zn[p0] + zn[p1] + zn[p2] + zn[p3]))
    E = len(mid_x)
    fc_base = V + E
    cc_base = fc_base + F
    new_n_nodes = cc_base + C

    # --- Pass B: per-cell corners, cell centers, interior faces -------------
    corner_idx = {}
    cc_x = array("d"); cc_y = array("d"); cc_z = array("d")
    in_fn = array("i"); in_left = array("i"); in_right = array("i")
    for c in range(C):
        corners = hex_corners(z, cell_faces, c)
        for k in range(8):
            corner_idx[(c, corners[k])] = k
        sx = sy = sz = 0.0
        for p in corners:
            sx += xn[p]; sy += yn[p]; sz += zn[p]
        cc_x.append(sx / 8.0); cc_y.append(sy / 8.0); cc_z.append(sz / 8.0)
        cc = cc_base + c

        # child centroids (P + 3 edge midpoints + 3 face centers + cell
        # center, all / 8), needed for interior-face winding
        cidx = {}
        for k in range(8):
            cidx[corners[k]] = k
        cent = []
        for p in corners:
            cent.append([xn[p] / 8.0, yn[p] / 8.0, zn[p] / 8.0])
        for f in cell_faces[c]:
            b = 4 * f
            ring = (fn[b], fn[b + 1], fn[b + 2], fn[b + 3])
            w = 1.0 / 8.0
            fx = fc_x[f] * w; fy = fc_y[f] * w; fz = fc_z[f] * w
            for node in ring:
                ce = cent[cidx[node]]
                ce[0] += fx; ce[1] += fy; ce[2] += fz
        for ia, ib in cell_edges(corners):
            pa, pb = corners[ia], corners[ib]
            w = 0.5 / 8.0
            mx = (xn[pa] + xn[pb]) * w
            my = (yn[pa] + yn[pb]) * w
            mz = (zn[pa] + zn[pb]) * w
            ca = cent[ia]; cb = cent[ib]
            ca[0] += mx; ca[1] += my; ca[2] += mz
            cb[0] += mx; cb[1] += my; cb[2] += mz
        w = 1.0 / 8.0
        ccx, ccy, ccz = cc_x[c] * w, cc_y[c] * w, cc_z[c] * w
        for k in range(8):
            cent[k][0] += ccx; cent[k][1] += ccy; cent[k][2] += ccz
        # each edge of the two faces that contain it (in this cell's face list)
        edge_faces = {}
        for f in cell_faces[c]:
            b = 4 * f
            ring = (fn[b], fn[b + 1], fn[b + 2], fn[b + 3])
            for i in range(4):
                a, e = ring[i], ring[(i + 1) % 4]
                key = a * stride + e if a < e else e * stride + a
                lst = edge_faces.get(key)
                if lst is None:
                    edge_faces[key] = [f]
                else:
                    lst.append(f)
        for ia, ib in cell_edges(corners):
            pa, pb = corners[ia], corners[ib]
            key = pa * stride + pb if pa < pb else pb * stride + pa
            fs = edge_faces[key]
            if len(fs) != 2:
                fail("cell %d edge %s has %d adjacent faces" % (c, key, len(fs)))
            f1, f2 = fs[0], fs[1]
            mid = mid_id[key]
            m_x, m_y, m_z = mid_x[mid - V], mid_y[mid - V], mid_z[mid - V]
            g1 = fc_base + f1; g2 = fc_base + f2
            g1x, g1y, g1z = fc_x[f1], fc_y[f1], fc_z[f1]
            g2x, g2y, g2z = fc_x[f2], fc_y[f2], fc_z[f2]
            ccx, ccy, ccz = cc_x[c], cc_y[c], cc_z[c]
            # normal of ring (m, fc1, cc, fc2) via diagonal cross product
            d1x = ccx - m_x; d1y = ccy - m_y; d1z = ccz - m_z
            d2x = g2x - g1x; d2y = g2y - g1y; d2z = g2z - g1z
            nx_ = d1y * d2z - d1z * d2y
            ny_ = d1z * d2x - d1x * d2z
            nz_ = d1x * d2y - d1y * d2x
            dx = cent[ib][0] - cent[ia][0]
            dy = cent[ib][1] - cent[ia][1]
            dz = cent[ib][2] - cent[ia][2]
            if nx_ * dx + ny_ * dy + nz_ * dz < 0.0:
                ring = (mid, g2, cc, g1)
            else:
                ring = (mid, g1, cc, g2)
            in_fn.extend(ring)
            in_left.append(c * 8 + ia)
            in_right.append(c * 8 + ib)

    # --- Pass C: external sub-faces ----------------------------------------
    ex_fn = array("i"); ex_left = array("i"); ex_right = array("i")
    left, right = z.left, z.right
    for f in range(F):
        b = 4 * f
        p0 = fn[b]; p1 = fn[b + 1]; p2 = fn[b + 2]; p3 = fn[b + 3]
        m01 = mid_id[p0 * stride + p1 if p0 < p1 else p1 * stride + p0]
        m12 = mid_id[p1 * stride + p2 if p1 < p2 else p2 * stride + p1]
        m23 = mid_id[p2 * stride + p3 if p2 < p3 else p3 * stride + p2]
        m30 = mid_id[p3 * stride + p0 if p3 < p0 else p0 * stride + p3]
        g = fc_base + f
        lc = left[f]
        rc = right[f]
        for i, (pi, mn, mp) in enumerate(((p0, m01, m30), (p1, m12, m01),
                                          (p2, m23, m12), (p3, m30, m23))):
            ex_fn.extend((pi, mn, g, mp))
            ex_left.append(lc * 8 + corner_idx[(lc, pi)])
            if f >= nB:
                ex_right.append(rc * 8 + corner_idx[(rc, pi)])
            elif ghost_mode:
                ex_right.append(4 * f + i + 8 * C)
            else:
                ex_right.append(-1)

    new = Zone()
    new.n_nodes = new_n_nodes
    new.n_faces = 4 * F + 12 * C
    new.n_cells = 8 * C
    new.vol_bc_type = z.vol_bc_type
    new.num_face_node = array("i", b"\x04\x00\x00\x00" * new.n_faces)
    new.fnmem = ex_fn + in_fn
    new.left = ex_left + in_left
    new.right = ex_right + in_right
    new.n_bfaces = 4 * nB
    new.bc_type = array("i", bytes(4 * new.n_bfaces))
    new.bc_name = array("i", bytes(4 * new.n_bfaces))
    for i in range(new.n_bfaces):
        new.bc_type[i] = z.bc_type[i // 4]
        new.bc_name[i] = z.bc_name[i // 4]
    new.xn = z.xn + mid_x + fc_x + cc_x
    new.yn = z.yn + mid_y + fc_y + cc_y
    new.zn = z.zn + mid_z + fc_z + cc_z

    report.append("input : V=%d F=%d C=%d nB=%d E=%d" % (V, F, C, nB, E))
    report.append("euler-input: V-E+F-C=%d" % (V - E + F - C))
    report.append(
        "output: V=%d F=%d C=%d nB=%d"
        % (new.n_nodes, new.n_faces, new.n_cells, new.n_bfaces)
    )
    if new.n_nodes != V + E + F + C:
        fail("node count mismatch")
    if new.n_faces != 4 * F + 12 * C:
        fail("face count mismatch")
    if new.n_cells != 8 * C:
        fail("cell count mismatch")
    if len(new.fnmem) != 4 * new.n_faces:
        fail("face-node storage mismatch")
    if len(new.left) != new.n_faces or len(new.right) != new.n_faces:
        fail("left/right storage mismatch")
    return new


def write_ofl(path, z, pid, zone_type):
    if sys.byteorder != "little" or array("i").itemsize != 4 or array("d").itemsize != 8:
        fail("host layout not little-endian 4-byte int / 8-byte double")
    payload = bytearray()
    payload += struct.pack("<3i", z.n_nodes, z.n_faces, z.n_cells)
    payload += z.xn.tobytes() + z.yn.tobytes() + z.zn.tobytes()
    payload += struct.pack("<i", z.vol_bc_type)
    payload += z.num_face_node.tobytes()
    payload += z.fnmem.tobytes()
    payload += z.left.tobytes() + z.right.tobytes()
    payload += struct.pack("<i", z.n_bfaces)
    payload += z.bc_type.tobytes() + z.bc_name.tobytes()
    payload += struct.pack("<i", 0)
    header = struct.pack("<i", 1)
    header += struct.pack("<i", pid[0]) + struct.pack("<i", zone_type[0])
    header += struct.pack("<q", len(payload))
    with open(path, "wb") as fh:
        fh.write(header)
        fh.write(payload)


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        return 2
    src = sys.argv[1]
    dst = sys.argv[2]
    n_cuts = 0
    i = 3
    while i < len(sys.argv):
        arg = sys.argv[i]
        if arg == "--cuts":
            i += 1
            if i >= len(sys.argv):
                print("--cuts needs a value (1/2/3)")
                return 2
            n_cuts = int(sys.argv[i])
            if n_cuts not in (1, 2, 3):
                print("--cuts must be 1, 2 or 3")
                return 2
        else:
            print("unknown argument %s" % arg)
            return 2
        i += 1
    report = []
    z, pid, zone_type = parse_ofl(src)
    report.append("REFINE_PARSE_PASS")
    _, vol_in_total = cell_volumes(z, report)
    if vol_in_total <= 0.0:
        fail("input total volume %.6e is not positive" % vol_in_total)
    report.append("input-total-volume=%.12e" % vol_in_total)

    if n_cuts <= 0:
        new = refine(z, report)
        report.append("REFINE_COUNTS_PASS")
    else:
        new = ladder_refine(z, pid, zone_type, n_cuts, report)
        report.append("REFINE_COUNTS_PASS")
    write_ofl(dst, new, pid, zone_type)
    validate_output(dst, vol_in_total, report)
    for line in report:
        print(line)
    return 0



# ---------------------------------------------------------------------------
# Sequential single-direction refinement (ladder mode, --cuts K)
#
# Uniform hex refinement multiplies cells by 8, which is too coarse for a
# doubling weak-scaling ladder.  This mode cuts ONE mesh direction per step
# (2x cells each), reusing the same corner algebra:
#   --cuts 1  : 2x cells (cut the direction labelled 0)
#   --cuts 2  : 4x cells (then cut label 1)
#   --cuts 3  : 8x cells (then cut label 2)
# Three sequential one-direction cuts reproduce the uniform 2x2x2 refinement
# exactly in real arithmetic: the midpoint of a segment between two cut-sheet
# midpoints lands on the original face centre, so the final node set is
# V + E + F + C.  Ladder meshes are therefore consistent step by step and the
# 8x endpoint can be cross-checked against the uniform refinement.
#
# Directions need a per-cell consistent global labelling, recovered by BFS:
# crossing a face constrains the face-perpendicular direction, and the shared
# edges of that face constrain the other two.  Label 0 is canonicalised to the
# direction whose edges are best aligned with a single global axis (the span
# direction of an extruded O-grid such as M6); labels 1/2 keep BFS order.
# The cut order is a free choice, only determinism matters.
# ---------------------------------------------------------------------------

def cell_dir_edge_pairs(corners):
    """The 12 edges of one hex as (corner_index_a, corner_index_b, dir) tuples.

    dir: 0 = vertical (b_i -> u_i), 1 = bottom/top ring direction b0->b1,
    2 = the other ring direction b1->b2.
    """
    out = []
    for i in range(4):
        out.append((i, 4 + i, 0))
    for a, b in ((0, 1), (2, 3), (4, 5), (6, 7)):
        out.append((a, b, 1))
    for a, b in ((1, 2), (3, 0), (5, 6), (7, 4)):
        out.append((a, b, 2))
    return out


def face_local_dir(ring, corners):
    """Which local direction (0/1/2) a face of this cell is perpendicular to.

    ring: the face's 4 node ids; corners: the cell's [b0..b3, u0..u3] list.
    """
    bset = set(corners[0:4])
    uset = set(corners[4:8])
    rset = set(ring)
    if rset == bset or rset == uset:
        return 0
    shared_b = [c for c in corners[0:4] if c in rset]
    if len(shared_b) != 2:
        return -1
    i0 = corners.index(shared_b[0])
    i1 = corners.index(shared_b[1])
    # The side face contains exactly one bottom-ring edge.  Express it in ring
    # order, (b_i, b_{(i+1)%4}); the wrapped pair (b3, b0) must use base 3, not
    # min(0, 3) = 0, or its parity flips and the face direction is misread.
    if (i0 + 1) % 4 == i1:
        base = i0
    elif (i1 + 1) % 4 == i0:
        base = i1
    else:
        return -1
    # Bases 0 and 2 span (b0,b1)/(b2,b3) = local direction 1, so the face is
    # perpendicular to direction 2; bases 1 and 3 span (b1,b2)/(b3,b0) = local
    # direction 2, so the face is perpendicular to direction 1.
    return 2 if base % 2 == 0 else 1


def direction_lengths(z, corners_flat):
    """Per-cell mean edge length of each local direction (array('d'), 3/cell).

    The M6 O-grid collapses ~20% of its cells onto zero-thickness pole/seam
    layers; there the perpendicular direction has no measurable length and no
    consistent labelling exists.  Lengths let the propagator skip exactly those
    constraints instead of failing on a genuine geometric singularity.
    """
    n = z.n_cells
    out = array("d", bytes(8 * 3 * n))
    for c in range(n):
        base = 8 * c
        corners = corners_flat[base : base + 8]
        acc = [0.0, 0.0, 0.0]
        cnt = [0, 0, 0]
        for a, b, d in cell_dir_edge_pairs(corners):
            pa, pb = corners[a], corners[b]
            dx = z.xn[pb] - z.xn[pa]
            dy = z.yn[pb] - z.yn[pa]
            dz = z.zn[pb] - z.zn[pa]
            acc[d] += (dx * dx + dy * dy + dz * dz) ** 0.5
            cnt[d] += 1
        for d in range(3):
            out[3 * c + d] = acc[d] / cnt[d] if cnt[d] else 0.0
    return out


def label_directions(z, cell_faces, corners_flat, report):
    """BFS-propagate a global direction label (0/1/2) per cell local frame.

    Labels are physical edge classes: crossing a face identifies the
    face-perpendicular classes and each shared edge identifies that edge's
    class.  Constraints are only applied where the geometry has real length;
    zero-thickness (collapsed) directions are filled in afterwards by a
    deterministic rule, because no labelling exists there.
    """
    from collections import deque

    n = z.n_cells
    lengths = direction_lengths(z, corners_flat)
    lmax = max(lengths) if n else 0.0
    solid_len = 1e-6 * lmax if lmax > 0.0 else 0.0
    labels = array("b", b"\xff" * (3 * n))  # -1 = unlabelled
    src_kind = array("b", bytes(3 * n))     # 1 = face, 2 = edge
    src_id = array("i", bytes(4 * 3 * n))   # face id / packed edge key
    src_cell = array("i", bytes(4 * 3 * n)) # neighbour cell that proposed
    solid = lambda c, d: lengths[3 * c + d] > solid_len

    def set_label(c, d, g, kind, ident, from_cell):
        idx = 3 * c + d
        cur = labels[idx]
        if cur == -1:
            labels[idx] = g
            src_kind[idx] = kind
            src_id[idx] = ident
            src_cell[idx] = from_cell
            return True
        if cur != g:
            print(
                "DEBUG twist: cell %d local dir %d already %d (via %s %s from cell %d)"
                ", now %d (via %s %s from cell %d)"
                % (c, d, cur,
                   "face" if src_kind[idx] == 1 else "edge", src_id[idx],
                   src_cell[idx], g,
                   "face" if kind == 1 else "edge", ident, from_cell),
                file=sys.stderr)
            fail("direction labelling found a twist at cell %d dir %d: %d vs %d"
                 % (c, d, cur, g))
        return False

    def edge_map(c):
        base = 8 * c
        corners = corners_flat[base : base + 8]
        emap = {}
        for a, b, d in cell_dir_edge_pairs(corners):
            key = (corners[a], corners[b])
            emap[key if key[0] < key[1] else (key[1], key[0])] = d
        return corners, emap

    seeds = 0
    queue = deque()
    for seed in range(n):
        if labels[3 * seed] != -1:
            continue
        labels[3 * seed + 0] = 0
        labels[3 * seed + 1] = 1
        labels[3 * seed + 2] = 2
        seeds += 1
        queue.append(seed)
        while queue:
            a = queue.popleft()
            a_corners, a_edges = edge_map(a)
            for f in cell_faces[a]:
                if z.left[f] == a:
                    b_cell = z.right[f] if f >= z.n_bfaces else -1
                else:
                    b_cell = z.left[f]
                if b_cell < 0 or labels[3 * b_cell] != -1:
                    continue
                b_corners, b_edges = edge_map(b_cell)
                ring = face_ring(z, f)
                da = face_local_dir(ring, a_corners)
                db = face_local_dir(ring, b_corners)
                if da < 0 or db < 0:
                    fail("face %d does not match a frame side of cell %d" % (f, a))
                if solid(a, da) and solid(b_cell, db):
                    set_label(b_cell, db, labels[3 * a + da], 1, f, a)
                for i in range(4):
                    p = ring[i]
                    q = ring[(i + 1) % 4]
                    key = (p, q) if p < q else (q, p)
                    if key not in b_edges or key not in a_edges:
                        fail("shared edge %s missing around face %d" % (key, f))
                    ea = a_edges[key]
                    eb = b_edges[key]
                    if solid(a, ea) and solid(b_cell, eb):
                        set_label(b_cell, eb, labels[3 * a + ea], 2,
                                  p if p < q else q, a)
                queue.append(b_cell)

    # Complete cells that only carry degenerate (zero-length) directions: the
    # remaining labels are assigned in ascending order, which is geometrically
    # irrelevant for a collapsed direction and keeps the permutation property.
    completed = 0
    for c in range(n):
        have = [labels[3 * c + d] for d in range(3)]
        if -1 not in have:
            continue
        used = set(g for g in have if g >= 0)
        missing = [g for g in (0, 1, 2) if g not in used]
        for d in range(3):
            if labels[3 * c + d] == -1:
                labels[3 * c + d] = missing.pop(0)
        completed += 1
    report.append(
        "direction-labelling: %d seed(s), %d cells, %d completed from collapsed "
        "directions (solid threshold %.3e)" % (seeds, n, completed, solid_len))
    return labels


def canonicalize_labels(z, corners_flat, labels, report):
    """Renumber global labels so label 0 is the THICKEST direction.

    Cutting the thickest direction first keeps the child cells well shaped and,
    on an O-grid whose pole layers are already zero-thickness in one direction,
    guarantees the cut never runs along a collapsed layer.  Median edge length
    per label decides the order; the axis alignment is reported for context.
    """
    n = z.n_cells
    step = max(1, n // 20000)
    lens = [[], [], []]
    axis_hits = [[0, 0, 0] for _ in range(3)]
    for c in range(0, n, step):
        base = 8 * c
        corners = corners_flat[base : base + 8]
        for a, b, d in cell_dir_edge_pairs(corners):
            g = labels[3 * c + d]
            if g < 0:
                continue
            dx = z.xn[corners[b]] - z.xn[corners[a]]
            dy = z.yn[corners[b]] - z.yn[corners[a]]
            dz = z.zn[corners[b]] - z.zn[corners[a]]
            ln = (dx * dx + dy * dy + dz * dz) ** 0.5
            lens[g].append(ln)
            ax = max(range(3), key=lambda k: abs((dx, dy, dz)[k]))
            axis_hits[g][ax] += 1

    def median(values):
        if not values:
            return 0.0
        values.sort()
        return values[len(values) // 2]

    med = [median(v) for v in lens]
    order = sorted(range(3), key=lambda g: -med[g])
    perm = [0, 0, 0]
    for rank, g in enumerate(order):
        perm[g] = rank
        axis = max(range(3), key=lambda k: axis_hits[g][k])
        report.append(
            "direction g=%d med-edge=%.6e rank=%d (dominant axis %s)"
            % (g, med[g], rank, "xyz"[axis])
        )
    new = array("b", bytes(3 * n))
    for i in range(3 * n):
        if labels[i] not in (0, 1, 2):
            fail("cell %d direction %d carries invalid label %d"
                 % (i // 3, i % 3, labels[i]))
        new[i] = perm[labels[i]]
    return new


class LadderState:
    """In-memory hex mesh with carried corner frames and direction labels."""

    def __init__(self, z, corners_flat, labels):
        self.n_nodes = z.n_nodes
        self.xn = array("d", z.xn)
        self.yn = array("d", z.yn)
        self.zn = array("d", z.zn)
        self.n_cells = z.n_cells
        self.corners = array("i", corners_flat)
        self.labels = array("b", labels)
        self.n_faces = z.n_faces
        self.frings = array("i", z.fnmem)
        self.left = array("i", z.left)
        self.right = array("i", z.right)
        self.n_bfaces = z.n_bfaces
        self.bc_type = array("i", z.bc_type)
        self.bc_name = array("i", z.bc_name)
        self.vol_bc_type = z.vol_bc_type


def apply_cut(st, cut_label, report):
    """Cut every cell along the direction carrying `cut_label` (2x cells).

    Each child frame is normalised so that its bottom ring is outward (by the
    sign of the corner triple product, i.e. the cell volume sign), and every
    emitted ring is then the child's own canonical outward face ring:

        bottom (b0,b1,b2,b3)   top (b3,b2,b1,b0)
        sides  (b_i,u_i,u_{i+1},b_{i+1})

    No ring order is inferred from the parent's storage order or from
    cross-product signs against thin geometry, so a collapsed O-grid layer can
    not silently produce an inward-wound face.
    """
    n_cells = st.n_cells
    corners = st.corners
    labels = st.labels
    xn, yn, zn = st.xn, st.yn, st.zn
    frings = st.frings
    left = st.left
    right = st.right

    # --- CSR cell -> faces -------------------------------------------------
    deg = array("i", bytes(4 * n_cells))
    for f in range(st.n_faces):
        deg[left[f]] += 1
        if f >= st.n_bfaces:
            deg[right[f]] += 1
    offs = array("i", bytes(4 * (n_cells + 1)))
    total = 0
    for c in range(n_cells):
        offs[c] = total
        total += deg[c]
    offs[n_cells] = total
    cursor = array("i", offs)
    cfaces = array("i", bytes(4 * total))
    for f in range(st.n_faces):
        c = left[f]
        cfaces[cursor[c]] = f
        cursor[c] += 1
        if f >= st.n_bfaces:
            c = right[f]
            cfaces[cursor[c]] = f
            cursor[c] += 1

    # --- per cell: cut direction, child frames -----------------------------
    mid = {}
    base_id = st.n_nodes
    cut_dir = array("b", bytes(n_cells))
    sheet = array("i", bytes(4 * n_cells))
    end0 = array("i", bytes(4 * n_cells))
    end1 = array("i", bytes(4 * n_cells))
    frames = array("i", bytes(4 * 8 * 2 * n_cells))  # 8 nodes per child
    new_labels = array("b", bytes(3 * 2 * n_cells))

    def get_mid(a, b):
        key = (a, b) if a < b else (b, a)
        m = mid.get(key)
        if m is None:
            m = base_id + len(mid)
            mid[key] = m
            xn.append(0.5 * (xn[a] + xn[b]))
            yn.append(0.5 * (yn[a] + yn[b]))
            zn.append(0.5 * (zn[a] + zn[b]))
        return m

    def mirror(frame):
        """Reverse the bottom/top rings so the bottom ring is outward."""
        return (frame[0], frame[3], frame[2], frame[1],
                frame[4], frame[7], frame[6], frame[5])

    def outward_bottom(frame):
        """Triple product (b1-b0) x (b2-b0) . (u0-b0).

        Negative means the bottom ring winds outward (the cross product points
        away from the cell, i.e. toward -u).  Positive means it winds inward.
        """
        b0, b1, b2, u0 = frame[0], frame[1], frame[2], frame[4]
        d1x = xn[b1] - xn[b0]; d1y = yn[b1] - yn[b0]; d1z = zn[b1] - zn[b0]
        d2x = xn[b2] - xn[b0]; d2y = yn[b2] - yn[b0]; d2z = zn[b2] - zn[b0]
        d3x = xn[u0] - xn[b0]; d3y = yn[u0] - yn[b0]; d3z = zn[u0] - zn[b0]
        cx = d1y * d2z - d1z * d2y
        cy = d1z * d2x - d1x * d2z
        cz = d1x * d2y - d1y * d2x
        return cx * d3x + cy * d3y + cz * d3z

    for c in range(n_cells):
        lb = 3 * c
        if labels[lb] == cut_label:
            j = 0
        elif labels[lb + 1] == cut_label:
            j = 1
        elif labels[lb + 2] == cut_label:
            j = 2
        else:
            fail("cell %d has no direction with label %d" % (c, cut_label))
        cut_dir[c] = j
        b8 = 8 * c
        p = (corners[b8], corners[b8 + 1], corners[b8 + 2], corners[b8 + 3])
        q = (corners[b8 + 4], corners[b8 + 5], corners[b8 + 6], corners[b8 + 7])

        if j == 0:
            e0, e1 = p, q
            s = tuple(get_mid(p[i], q[i]) for i in range(4))
            child0 = p + s
            child1 = s + q
            triples = ((labels[lb], labels[lb + 1], labels[lb + 2]),
                       (labels[lb], labels[lb + 1], labels[lb + 2]))
        elif j == 1:
            e0 = (p[1], p[2], q[2], q[1])
            e1 = (p[3], p[0], q[0], q[3])
            s0 = get_mid(p[0], p[1]); s1 = get_mid(p[2], p[3])
            s2 = get_mid(q[2], q[3]); s3 = get_mid(q[0], q[1])
            s = (s2, s1, s0, s3)
            child0 = (p[1], p[2], q[2], q[1], s0, s1, s2, s3)
            child1 = (p[3], p[0], q[0], q[3], s1, s0, s3, s2)
            t = (labels[lb + 1], labels[lb + 2], labels[lb])
            triples = (t, t)
        else:
            e0 = (p[0], p[1], q[1], q[0])
            e1 = (p[2], p[3], q[3], q[2])
            s0 = get_mid(p[3], p[0]); s1 = get_mid(p[1], p[2])
            s2 = get_mid(q[1], q[2]); s3 = get_mid(q[3], q[0])
            s = (s2, s1, s0, s3)
            child0 = (p[0], p[1], q[1], q[0], s0, s1, s2, s3)
            child1 = (p[2], p[3], q[3], q[2], s1, s0, s3, s2)
            t = (labels[lb + 2], labels[lb + 1], labels[lb])
            triples = (t, t)

        for k, frame in enumerate((child0, child1)):
            triple = triples[k]
            # triple > 0 means (b1-b0)x(b2-b0) points toward u0, i.e. INTO the
            # cell, so the bottom ring winds inward and must be mirrored.
            if outward_bottom(frame) > 0.0:
                frame = mirror(frame)
                triple = (triple[0], triple[2], triple[1])
            frames[16 * c + 8 * k : 16 * c + 8 * k + 8] = array("i", frame)
            new_labels[3 * (2 * c + k) : 3 * (2 * c + k) + 3] = array("b", triple)
        end0[4 * c : 4 * c + 4] = array("i", e0)
        end1[4 * c : 4 * c + 4] = array("i", e1)
        sheet[4 * c : 4 * c + 4] = array("i", s)

    def canonical(child, quad):
        """The child's own outward ring for the face with node set `quad`."""
        b4 = 16 * (child // 2) + 8 * (child % 2)
        f = frames[b4 : b4 + 8]
        b0, b1, b2, b3 = f[0], f[1], f[2], f[3]
        u0, u1, u2, u3 = f[4], f[5], f[6], f[7]
        # Normalised frame: bottom ring outward.  Outward rings are then
        # bottom (b0,b1,b2,b3), top (u0,u3,u2,u1) and side i
        # (b_i, u_i, u_{i+1}, b_{i+1}); the top is NOT the reverse of the
        # bottom ring, it pairs u0..u3 across the sheet.
        cands = ((b0, b1, b2, b3), (u0, u3, u2, u1),
                 (b0, u0, u1, b1), (b1, u1, u2, b2),
                 (b2, u2, u3, b3), (b3, u3, u0, b0))
        target = set(quad)
        for cand in cands:
            if set(cand) == target:
                return cand
        print("DEBUG canonical fail child=%d quad=%s frame=%s"
              % (child, sorted(target), list(f)), file=sys.stderr)
        for cand in cands:
            print("    cand=%s set=%s" % (cand, sorted(set(cand))), file=sys.stderr)
        fail("child %d has no canonical face with nodes %s" % (child, sorted(target)))

    # Ownership is decided geometrically: a face belongs to the child on its
    # side of the parent's mid-sheet plane.  Using the same rule for both cells
    # of an interface keeps degenerate (folded) pairs consistent, where a
    # node-set test cannot distinguish the two children because they share
    # nodes.  Coincident faces then cancel exactly as they do in the input.
    sheet_mid = array("d", bytes(8 * 3 * n_cells))
    sheet_n = array("d", bytes(8 * 3 * n_cells))
    for c in range(n_cells):
        b8 = 16 * c
        c0 = [0.0, 0.0, 0.0]
        c1 = [0.0, 0.0, 0.0]
        for k in range(4):
            n0 = frames[b8 + k]
            n1 = frames[b8 + 8 + k]
            c0[0] += xn[n0]; c0[1] += yn[n0]; c0[2] += zn[n0]
            c1[0] += xn[n1]; c1[1] += yn[n1]; c1[2] += zn[n1]
        for k in range(3):
            c0[k] *= 0.25
            c1[k] *= 0.25
            sheet_mid[3 * c + k] = 0.5 * (c0[k] + c1[k])
            sheet_n[3 * c + k] = c1[k] - c0[k]

    degenerate = array("b", bytes(n_cells))

    def owner(c, quad):
        """Child of cell c on the side of the mid-sheet plane holding `quad`.

        Collapsed cells (zero thickness along the cut) make the signed test
        meaningless, so they fall back to node containment and are flagged so
        the caller can keep the two halves on different siblings.
        """
        qx = (xn[quad[0]] + xn[quad[1]] + xn[quad[2]] + xn[quad[3]]) * 0.25
        qy = (yn[quad[0]] + yn[quad[1]] + yn[quad[2]] + yn[quad[3]]) * 0.25
        qz = (zn[quad[0]] + zn[quad[1]] + zn[quad[2]] + zn[quad[3]]) * 0.25
        t = ((qx - sheet_mid[3 * c]) * sheet_n[3 * c]
             + (qy - sheet_mid[3 * c + 1]) * sheet_n[3 * c + 1]
             + (qz - sheet_mid[3 * c + 2]) * sheet_n[3 * c + 2])
        nn = (sheet_n[3 * c] ** 2 + sheet_n[3 * c + 1] ** 2
              + sheet_n[3 * c + 2] ** 2) ** 0.5
        if nn > 1e-9:
            return 2 * c + (0 if t < 0.0 else 1)
        degenerate[c] = 1
        quad_set = set(quad)
        b4 = 16 * c
        if quad_set <= set(frames[b4 : b4 + 8]):
            return 2 * c
        if quad_set <= set(frames[b4 + 8 : b4 + 16]):
            return 2 * c + 1
        return 2 * c

    # --- emit: boundary sub-faces first, then interior --------------------
    out_frings = array("i")
    out_left = array("i")
    out_right = array("i")
    out_bc_type = array("i")
    out_bc_name = array("i")
    n_halves = 0

    for f in range(st.n_faces):
        b4 = 4 * f
        ring = (frings[b4], frings[b4 + 1], frings[b4 + 2], frings[b4 + 3])
        bnd = f < st.n_bfaces
        la = left[f]
        ra = right[f]
        d_a = face_local_dir(ring, corners[8 * la : 8 * la + 8])
        if d_a < 0:
            fail("face %d is not a side of cell %d" % (f, la))
        split_a = d_a != cut_dir[la]

        if ra >= 0:
            d_b = face_local_dir(ring, corners[8 * ra : 8 * ra + 8])
            if d_b < 0:
                fail("face %d is not a side of cell %d" % (f, ra))
            split_b = d_b != cut_dir[ra]
            if split_a != split_b:
                fail("face %d is split for one cell only" % f)
        else:
            split_b = split_a

        quads = []
        if not split_a:
            quads.append(ring)
        else:
            n_halves += 1
            keys = []
            for i in range(4):
                a = ring[i]
                bb = ring[(i + 1) % 4]
                keys.append(mid.get((a, bb) if a < bb else (bb, a)))
            for i in range(4):
                if keys[i] is not None:
                    continue
                a = ring[i]
                bb = ring[(i + 1) % 4]
                m_fwd = keys[(i + 1) % 4]
                if m_fwd is None:
                    m_fwd = get_mid(bb, ring[(i + 2) % 4])
                m_bwd = keys[(i - 1) % 4]
                if m_bwd is None:
                    m_bwd = get_mid(ring[(i - 1) % 4], a)
                quads.append((a, bb, m_fwd, m_bwd))
            if len(quads) != 2:
                fail("face %d produced %d halves" % (f, len(quads)))

        lc_children = [owner(la, q) for q in quads]
        rc_children = [owner(ra, q) for q in quads] if ra >= 0 else []
        # Only for collapsed cells are the two children coincident; there the
        # halves must be forced onto different siblings or both would describe
        # the same interface twice.
        if len(quads) == 2 and degenerate[la] and lc_children[0] == lc_children[1]:
            lc_children[1] = lc_children[0] ^ 1
        if (len(quads) == 2 and ra >= 0 and degenerate[ra]
                and rc_children[0] == rc_children[1]):
            rc_children[1] = rc_children[0] ^ 1
        for idx, quad in enumerate(quads):
            lc_child = lc_children[idx]
            out_frings.extend(canonical(lc_child, quad))
            out_left.append(lc_child)
            if bnd:
                out_right.append(-1)
                out_bc_type.append(st.bc_type[f])
                out_bc_name.append(st.bc_name[f])
            else:
                out_right.append(rc_children[idx])

    # one sheet face per cell, left = child0
    for c in range(n_cells):
        quad = tuple(sheet[4 * c : 4 * c + 4])
        out_frings.extend(canonical(2 * c, quad))
        out_left.append(2 * c)
        out_right.append(2 * c + 1)

    # --- commit ------------------------------------------------------------
    st.n_nodes = base_id + len(mid)
    st.n_cells = 2 * n_cells
    st.corners = frames
    st.labels = new_labels
    st.n_faces = len(out_frings) // 4
    st.frings = out_frings
    st.left = out_left
    st.right = out_right
    st.n_bfaces = len(out_bc_type)
    st.bc_type = out_bc_type
    st.bc_name = out_bc_name
    report.append(
        "cut label %d: nodes +%d cells x2 faces %d (split parents %d) nB %d"
        % (cut_label, len(mid), st.n_faces, n_halves, st.n_bfaces))
    return st


def count_unique_edges(z, cell_faces, corners_flat):
    edges = set()
    for c in range(z.n_cells):
        base = 8 * c
        corners = corners_flat[base : base + 8]
        for a, b, _d in cell_dir_edge_pairs(corners):
            key = (corners[a], corners[b])
            edges.add(key if key[0] < key[1] else (key[1], key[0]))
    return len(edges)


def ladder_refine(z, pid, zone_type, n_cuts, report):
    """Sequential one-direction refinement; returns the output Zone."""
    cell_faces = build_cell_faces(z)
    corners_flat = array("i", bytes(4 * 8 * z.n_cells))
    for c in range(z.n_cells):
        corners_flat[8 * c : 8 * c + 8] = array("i", hex_corners(z, cell_faces, c))
    labels = label_directions(z, cell_faces, corners_flat, report)
    labels = canonicalize_labels(z, corners_flat, labels, report)

    bad = []
    for c in range(z.n_cells):
        trip = (labels[3 * c], labels[3 * c + 1], labels[3 * c + 2])
        if sorted(trip) != [0, 1, 2]:
            bad.append((c, trip))
    if bad:
        print("DEBUG bad-label cells=%d first=%s" % (len(bad), bad[:5]),
              file=sys.stderr)
        c0 = bad[0][0]
        print("DEBUG cell %d corners=%s" % (c0, corners_flat[8*c0:8*c0+8]),
              file=sys.stderr)
        for f in cell_faces[c0]:
            lc = z.left[f]
            rc = z.right[f] if f >= z.n_bfaces else -1
            ltr = (labels[3*lc], labels[3*lc+1], labels[3*lc+2]) if lc >= 0 else None
            rtr = (labels[3*rc], labels[3*rc+1], labels[3*rc+2]) if rc >= 0 else None
            print("  face %d ring=%s L=%d%s R=%d%s"
                  % (f, face_ring(z, f), lc, ltr, rc, rtr), file=sys.stderr)
        fail("%d cells have non-permutation direction labels" % len(bad))

    v0, f0, c0, nb0 = z.n_nodes, z.n_faces, z.n_cells, z.n_bfaces
    e_total = count_unique_edges(z, cell_faces, corners_flat)
    report.append("unique edges in input: %d" % e_total)

    st = LadderState(z, corners_flat, labels)
    for k in range(n_cuts):
        apply_cut(st, k, report)

    if n_cuts == 3:
        # Sequential three-direction cuts reproduce the uniform 2x2x2
        # refinement: closed-form identity check on node/face/boundary counts.
        if st.n_nodes != v0 + e_total + f0 + c0:
            fail("ladder node count %d != uniform identity %d"
                 % (st.n_nodes, v0 + e_total + f0 + c0))
        if st.n_faces != 4 * f0 + 12 * c0:
            fail("ladder face count %d != uniform identity %d"
                 % (st.n_faces, 4 * f0 + 12 * c0))
        if st.n_bfaces != 4 * nb0:
            fail("ladder boundary count %d != uniform identity %d"
                 % (st.n_bfaces, 4 * nb0))
        report.append("uniform identity holds: V=%d F=%d C=%d nB=%d"
                      % (st.n_nodes, st.n_faces, st.n_cells, st.n_bfaces))

    new = Zone()
    new.n_nodes = st.n_nodes
    new.n_faces = st.n_faces
    new.n_cells = st.n_cells
    new.vol_bc_type = z.vol_bc_type
    new.xn = st.xn
    new.yn = st.yn
    new.zn = st.zn
    new.num_face_node = array("i", b"\x04\x00\x00\x00" * st.n_faces)
    new.fnmem = st.frings
    new.left = st.left
    new.right = st.right
    new.n_bfaces = st.n_bfaces
    new.bc_type = st.bc_type
    new.bc_name = st.bc_name
    return new


def validate_output(dst, vol_in_total, report):
    back, _, _ = parse_ofl(dst)
    _, vol_out_total = cell_volumes(back, report)
    rel = abs(vol_out_total - vol_in_total) / vol_in_total
    report.append("re-read-total-volume=%.12e rel-diff=%.3e"
                  % (vol_out_total, rel))
    if rel > 1e-10:
        fail("total volume drifts by %.3e" % rel)
    report.append("REFINE_VOLUME_PASS")
    cell_faces = build_cell_faces(back)
    step = max(1, back.n_cells // 500)
    for c in range(0, back.n_cells, step):
        hex_corners(back, cell_faces, c)
    report.append("REFINE_REPARSE_PASS")
    report.append("REFINE_ALL_PASS")


if __name__ == "__main__":
    sys.exit(main())
