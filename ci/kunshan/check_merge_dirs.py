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
    # Node-identity correspondence: the SAME node id n must satisfy
    # latL[n] == T(latR[n]) + offset on all four ring nodes.  The ring start
    # may rotate, but a node may never be remapped to a different node whose
    # coordinates happen to match.  T = permutation x signs (48 candidates);
    # the offset follows from any one node and is then verified on the rest.
    import itertools
    for perm in itertools.permutations(range(3)):
        for signs in itertools.product((1, -1), repeat=3):
            rmap = {r: (perm[r], signs[r]) for r in range(3)}
            # offset from the first ring node, then verify all four BY ID
            n0 = ring[0]
            t0 = apply_rmap(latR[n0], rmap)
            off = tuple(latL[n0][j] - t0[j] for j in range(3))
            ok = True
            for n in ring:
                t = apply_rmap(latR[n], rmap)
                if tuple(t[j] + off[j] for j in range(3)) != latL[n]:
                    ok = False
                    break
            if ok:
                return (nL, nR, rmap), None
    return None, "no node-identity affine correspondence"


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
    """Deterministic exact solver for a bijection-constraint network.

    For each connected component, tries each root direction 0/1/2; BFS along
    forward and inverse permutations propagates labels edge by edge; any
    contradiction aborts that seed.  A component is UNSAT only when all three
    seeds fail, with the blocking path recorded per seed.  At most three graph
    traversals per component, no data structure beyond adjacency lists.
    """
    from collections import deque

    # verify each constraint is a bijection and extract P
    edges = []
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
        edges.append((f, l, r, tuple(P)))

    adj = [[] for _ in range(n_cells)]
    for (f, l, r, P) in edges:
        adj[l].append((f, r, P))
        adj[r].append((f, l, inverse(P)))

    seen = [False] * n_cells
    assign = [-1] * n_cells
    for root in range(n_cells):
        if seen[root]:
            continue
        comp = [root]
        seen[root] = True
        dq = deque([root])
        while dq:
            c = dq.popleft()
            for (f, o, P) in adj[c]:
                if not seen[o]:
                    seen[o] = True
                    comp.append(o)
                    dq.append(o)

        solved = None
        failures = []
        for seed in range(3):
            val = {root: seed}
            path = [root]
            block = None
            dq = deque([root])
            while dq and block is None:
                c = dq.popleft()
                for (f, o, P) in adj[c]:
                    if o == c:
                        # self-loop: consistency check
                        if P[val[c]] != val[c]:
                            block = (f, c, c, P)
                            break
                        continue
                    nv = P[val[c]]
                    if o in val:
                        if val[o] != nv:
                            block = (f, c, o, P)
                            break
                    else:
                        val[o] = nv
                        path.append(o)
                        dq.append(o)
            if block is None:
                # full re-check over every edge of the component
                for (f, l, r, P) in edges:
                    if l in val and r in val:
                        if P[val[l]] != val[r]:
                            block = (f, l, r, P)
                            break
            if block is None:
                solved = val
                break
            failures.append((seed, block))

        if solved is None:
            return ("UNSAT",
                    "component rooted at cell %d: all 3 root directions fail; "
                    "per-seed blocking edges: %s"
                    % (root, failures))
        for c, v in solved.items():
            assign[c] = v

    # final independent pass: every constraint verified against the assignment
    for (f, l, r, am) in cons:
        if not am[assign[l]][assign[r]]:
            return ("INCOMPLETE",
                    "post-check violation at face %d (%d,%d)" % (f, l, r))
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
          % (mode, dict(Counter(assign))))
    return 0


if __name__ == "__main__":
    sys.exit(main())
