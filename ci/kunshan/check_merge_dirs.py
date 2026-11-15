"""Candidate-direction compatibility for a per-parent merge (writes no grid).

Cross-parent octant differences are NOT a general necessity: the two parents of
a shared face may carry different local axis conventions.  This tool therefore
works only from real face identity and the 4-corner correspondence of every
shared parent face, and builds an explicit constraint system:

  variable per parent cell   : the split direction d in that cell's own frame
                               (2N: bisect along d -> 2 bricks of 4 fine cells;
                                4N: leave d unsplit -> 4 bricks of 2 cells)
  constraint per shared face : let aL/aR be the face-perpendicular axes in the
                               left/right frames and pi the in-face axis map
                               recovered from the ring corner correspondence.
                                 dL==aL and dR==aR            -> allowed
                                 exactly one equals its a     -> forbidden
                                 otherwise                    -> require pi(dL)==dR

Faces are never merged by unordered(left,right): two distinct parent faces can
join the same cell pair (the 72 O-grid pinch pairs) and each keeps its own
identity.  The solver propagates the domains; if a domain empties it reports the
offending face.  Verdicts are invariant under node/face relabelling because only
face identity and corner correspondence are used.

Usage: check_merge_dirs.py <parent.ofl> [--mode 2N|4N]
"""
import sys
from array import array
sys.path.insert(0, "ci/kunshan")
import refine_ofl as R
import agglomerate_ofl as A

CH = A.CHILD_LATTICE  # corner index -> (a,b,c) in the cell's own frame


DIR_TO_LAT = {0: 2, 1: 0, 2: 1}    # face_local_dir axis -> lattice (a,b,c)
LAT_TO_DIR = {2: 0, 0: 1, 1: 2}
CH = A.CHILD_LATTICE


def cell_lattice_positions(corners):
    """Lattice position of every corner node id in the cell's own frame.

    corners: the [b0..b3, u0..u3] list recovered by hex_corners.  The frame
    convention is fixed (b0->b1 is lattice a, b1->b2 is lattice b, b->u is
    lattice c), so corner k of index i maps to CHILD_LATTICE[i] -- by
    construction of refine_ofl, not by assumption about ring starts.
    """
    return {corners[i]: tuple(CH[i]) for i in range(8)}


def analyse_face(ring, latL, latR):
    """Full four-corner affine correspondence for one shared parent face.

    Returns (nL, nR, rmap) in LATTICE axes, or (None, reason).
    nL/nR: the face-perpendicular lattice axis on each side.  rmap maps a
    right lattice axis to (left axis, sign).  The offset is irrelevant for
    direction compatibility and is only used to verify the correspondence.
    """
    # node -> lattice coords per side
    ringL = [latL[n] for n in ring]
    ringR = [latR[n] for n in ring]
    if any(n not in latL for n in ring) or any(n not in latR for n in ring):
        return None, "ring node outside a cell frame"
    nL = lattice_normal(ringL)
    nR = lattice_normal(ringR)
    if nL is None or nR is None:
        return None, "face not lattice-planar"
    # enumerate axis maps: 3! permutations x 2^3 sign patterns
    import itertools
    for perm in itertools.permutations(range(3)):
        for signs in itertools.product((1, -1), repeat=3):
            rmap = {r: (perm[r], signs[r]) for r in range(3)}
            if affine_ok(ringR, ringL, rmap):
                return (nL, nR, rmap), None
    return None, "no affine corner correspondence"


def lattice_normal(ring):
    """Face-perpendicular lattice axis: edges span exactly two axes."""
    axes = set()
    for i in range(4):
        p, q = ring[i], ring[(i + 1) % 4]
        for j in range(3):
            if p[j] != q[j]:
                axes.add(j)
    if len(axes) != 2:
        return None
    return ({0, 1, 2} - axes).pop()


def affine_ok(ringR, ringL, rmap):
    """True if some offset makes rmap map ringR exactly onto ringL's set."""
    posL = set(ringL)
    for start in range(4):
        off = [ringL[0][j] - apply_rmap(ringR[start], rmap)[j] for j in range(3)]
        got = set()
        ok = True
        for r in ringR:
            m = apply_rmap(r, rmap)
            q = tuple(m[j] + off[j] for j in range(3))
            if q in posL:
                got.add(q)
            else:
                ok = False
                break
        if ok and len(got) == 4:
            return True
    return False


def apply_rmap(p, rmap):
    out = [0, 0, 0]
    for raxis in range(3):
        laxis, sign = rmap[raxis]
        out[laxis] = sign * p[raxis]
    return tuple(out)


def allowed_matrix(nL, nR, rmap, mode):
    """3x3 allowed combinations of lattice split axes (lattice indexing).

    2N: split along d -> bricks bisect along d.  A side cuts this face iff
        d != n.  If both cut, the cut lines coincide iff rmap maps dR to dL.
    4N: d is the UNSPLIT axis.  A side cuts this face iff d == n is FALSE...
        if d == n the face is perpendicular to the unsplit axis, so both
        in-face axes are split and the face is cut 2x2.  Compatibility: the
        set of cut axes must correspond through rmap.
    """
    allowed = [[0] * 3 for _ in range(3)]
    for dL in range(3):
        for dR in range(3):
            if mode == "2N":
                sL = dL != nL
                sR = dR != nR
                if sL != sR:
                    allowed[dL][dR] = 0
                elif not sL:
                    allowed[dL][dR] = 1
                else:
                    laxis, sign = rmap[dR]
                    allowed[dL][dR] = 1 if laxis == dL else 0
            else:
                # 4N: merge pairs of fine cells along d -> 4 bricks of 2.
                # Parent face with normal n: each brick contains children
                # adjacent to the face.  If d == n every brick holds exactly
                # one adjacent child, so the face splits into 4 sub-faces
                # (one per brick, fine-face granular).  If d != n every brick
                # holds two adjacent children, so the face splits into 2
                # pieces and the dividing lines coincide iff rmap(dR)==dL.
                # A whole (uncut) face never occurs in 4N.
                if dL == nL and dR == nR:
                    allowed[dL][dR] = 1
                elif (dL == nL) != (dR == nR):
                    allowed[dL][dR] = 0     # 4 vs 2 pieces: non-conformal
                else:
                    allowed[dL][dR] = 1 if rmap[dR][0] == dL else 0
    return allowed


def compose(p, q):
    return tuple(p[q[x]] for x in range(3))


def inverse(p):
    out = [0, 0, 0]
    for x in range(3):
        out[p[x]] = x
    return tuple(out)


def fixed_points(p):
    return frozenset(x for x in range(3) if p[x] == x)


def solve_perm_network(n_cells, cons):
    """Exact solver when every face constraint is a bijection dR = P(dL).

    Coloured union-find.  Invariant: perm[c] maps label(find(c)) -> label(c);
    path compression composes as perm[y] = perm[y] o perm[parent] (parent's
    mapping applied FIRST), i.e. label(y) = perm[y](label(root)).

    Constraint label(b) = p(label(a)) with a rooted at ra, b rooted at rb:
        label(b) = perm[b](label(rb)), label(a) = perm[a](label(ra))
        => perm[b](label(rb)) = p(perm[a](label(ra)))
    Attaching rb under ra needs perm[rb]_new: label(ra) -> label(rb):
        perm[rb]_new = inverse(perm[b]) o p o perm[a]
    Cycle edge (ra == rb) forces l = label(ra) to satisfy
        cyc(l) = l  with cyc = inverse(perm[b]) o p o perm[a]
    so the admissible root labels are the intersection of the fixed-point
    sets of all cycle permutations; empty intersection -> UNSAT with that
    cycle as witness.  Deterministic, near O(n alpha(n)); no search.
    """
    parent = list(range(n_cells))
    perm = [tuple(range(3)) for _ in range(n_cells)]
    size = [1] * n_cells
    cycles = []

    def find(x):
        root = x
        while parent[root] != root:
            root = parent[root]
        # first pass: compose perms bottom-up correctly
        while parent[x] != root:
            nxt = parent[x]
            perm[x] = compose(perm[x], perm[nxt])
            parent[x] = root
            x = nxt
        return root

    def union(a, b, p):
        ra, rb = find(a), find(b)
        pa, pb = perm[a], perm[b]
        if ra == rb:
            cyc = compose(inverse(pb), compose(p, pa))
            if cyc != tuple(range(3)):
                cycles.append((a, b, cyc))
            return
        shift = compose(inverse(pb), compose(p, pa))   # label(ra)->label(rb)
        if size[ra] < size[rb]:
            # attach ra under rb instead: need label(rb)->label(ra)
            parent[ra] = rb
            # perm[ra]_new: label(rb) -> label(ra); perm[ra]_old maps ra->ra
            perm[ra] = compose(inverse(shift), perm[ra])
            size[rb] += size[ra]
        else:
            parent[rb] = ra
            perm[rb] = compose(shift, perm[rb])
            size[ra] += size[rb]

    for (f, l, r, am) in cons:
        P = [None] * 3
        ok = True
        for dL in range(3):
            row = [dR for dR in range(3) if am[dL][dR]]
            if len(row) != 1:
                ok = False
                break
            P[dL] = row[0]
        if not ok or len(set(P)) != 3:
            return ("INCOMPLETE",
                    "face %d allowed matrix is not a bijection: %s" % (f, am))
        union(l, r, tuple(P))

    allowed = {}
    for c in range(n_cells):
        find(c)
        rt = parent[c]
        if rt not in allowed:
            allowed[rt] = frozenset(range(3))
    for (a, b, cyc) in cycles:
        rt = parent[a]
        allowed[rt] = allowed[rt] & fixed_points(cyc)
        if not allowed[rt]:
            return ("UNSAT",
                    "cycle (cells %d,%d) composes to permutation %s with no "
                    "fixed point: every candidate split direction is "
                    "excluded around this loop" % (a, b, cyc))
    assign = {}
    for c in range(n_cells):
        rt = parent[c]
        if not allowed.get(rt):
            return ("UNSAT", "component root %d has empty label set" % rt)
        assign[c] = perm[c][min(allowed[rt])]
    return ("PASS", assign)


def main():
    path = None
    mode = "2N"
    i = 1
    while i < len(sys.argv):
        if sys.argv[i] == "--mode":
            mode = sys.argv[i + 1]; i += 2
        else:
            path = sys.argv[i]; i += 1
    if path is None:
        print("usage: check_merge_dirs.py <parent.ofl> [--mode 2N|4N]"); return 2
    z, _, _ = R.parse_ofl(path)
    cf = R.build_cell_faces(z)
    lat = []
    for c in range(z.n_cells):
        corners = R.hex_corners(z, cf, c)
        lat.append(cell_lattice_positions(corners))

    cons = []
    invalid = []
    for f in range(z.n_bfaces, z.n_faces):
        l, r = z.left[f], z.right[f]
        got, err = analyse_face(R.face_ring(z, f), lat[l], lat[r])
        if got is None:
            invalid.append((f, err))
            continue
        nL, nR, rmap = got
        cons.append((f, l, r, allowed_matrix(nL, nR, rmap, mode)))
    print("mode=%s cells=%d interior faces=%d parsed=%d invalid=%d"
          % (mode, z.n_cells, z.n_faces - z.n_bfaces, len(cons), len(invalid)))
    if invalid:
        print("VERDICT=INVALID_INPUT")
        for f, err in invalid[:5]:
            print("  face %d: %s" % (f, err))
        return 2

    dom = [0b111] * z.n_cells
    by = [[] for _ in range(z.n_cells)]
    for (f, l, r, am) in cons:
        by[l].append((f, r, am, True))
        by[r].append((f, l, am, False))

    verdict, detail = solve_perm_network(z.n_cells, cons)
    if verdict == "INCOMPLETE":
        print("VERDICT=INCOMPLETE (%s)" % detail)
        return 3
    if verdict == "UNSAT":
        print("VERDICT=UNSAT")
        print("witness: %s" % detail)
        return 1
    assign = detail
    viol = 0
    for (f, l, r, am) in cons:
        if not am[assign[l]][assign[r]]:
            viol += 1
            if viol <= 3:
                print("  violation face %d pair (%d,%d)" % (f, l, r))
    if viol:
        print("VERDICT=INCOMPLETE (%d violations after assignment)" % viol)
        return 3
    from collections import Counter
    print("VERDICT=PASS mode=%s histogram=%s"
          % (mode, dict(Counter(assign.values()))))
    return 0


if __name__ == "__main__":
    sys.exit(main())
