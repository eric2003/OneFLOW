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


def face_axis_and_map(ring, cl, cr):
    """Perpendicular axis in each frame + in-face axis map (left -> right)."""
    aL = R.face_local_dir(ring, cl)
    aR = R.face_local_dir(ring, cr)
    if aL < 0 or aR < 0:
        return None
    oL = [x for x in range(3) if x != aL]
    oR = [x for x in range(3) if x != aR]
    # in-face coordinate of each ring node, per side
    def coords(cl, axes, node):
        k = cl.index(node)
        c = CH[k]
        return (c[axes[0]], c[axes[1]])
    # find the node at (1,0) on the left to learn the axis correspondence
    node10 = None
    for n in ring:
        if coords(cl, oL, n) == (1, 0):
            node10 = n
            break
    if node10 is None:
        return None
    uv = coords(cr, oR, node10)
    if uv == (1, 0):
        return aL, aR, oL, oR, oR[0], oR[1]  # uL->uR, vL->vR
    if uv == (0, 1):
        return aL, aR, oL, oR, oR[1], oR[0]  # uL->vR, vL->uR
    return None


def main():
    path = sys.argv[1]
    mode = "2N"
    if "--mode" in sys.argv:
        mode = sys.argv[sys.argv.index("--mode") + 1]
    z, _, _ = R.parse_ofl(path)
    cf = R.build_cell_faces(z)
    corners = [R.hex_corners(z, cf, c) for c in range(z.n_cells)]

    cons = []          # (left, right, allowed[3][3])
    bad_faces = []
    for f in range(z.n_bfaces, z.n_faces):
        l, r = z.left[f], z.right[f]
        ring = R.face_ring(z, f)
        got = face_axis_and_map(ring, corners[l], corners[r])
        if got is None:
            bad_faces.append(f)
            continue
        aL, aR, oL, oR, mL, mR = got
        # map left in-face axis -> right in-face axis
        axis_to_right = {}
        axis_to_right[oL[0]] = oR[0] if (mL == oR[0] and mR == oR[1]) else oR[1]
        axis_to_right[oL[1]] = oR[1] if (mL == oR[0] and mR == oR[1]) else oR[0]
        allowed = [[0] * 3 for _ in range(3)]
        for dL in range(3):
            for dR in range(3):
                if mode == "2N":
                    sL = (dL != aL)
                    sR = (dR != aR)
                    if not sL and not sR:
                        ok = True
                    elif sL != sR:
                        ok = False
                    else:
                        ok = axis_to_right[dL] == dR
                else:  # 4N: d is the UNSPLIT axis; the other two are split
                    sL = (dL == aL)
                    sR = (dR == aR)
                    if sL and sR:
                        ok = True
                    elif sL != sR:
                        ok = False
                    else:
                        othersL = [x for x in range(3) if x != dL]
                        othersR = [x for x in range(3) if x != dR]
                        ok = set(axis_to_right[x] for x in othersL) == set(othersR)
                allowed[dL][dR] = 1 if ok else 0
        if not any(any(row) for row in allowed):
            bad_faces.append(f)
        cons.append((l, r, allowed))

    print("mode=%s cells=%d shared faces=%d unanalysable faces=%d"
          % (mode, z.n_cells, len(cons), len(bad_faces)))
    if bad_faces:
        print("UNANALYSABLE first faces: %s" % bad_faces[:5])

    dom = [0b111] * z.n_cells
    by_cell = [[] for _ in range(z.n_cells)]
    for idx, (l, r, allowed) in enumerate(cons):
        by_cell[l].append((idx, r, allowed, True))
        by_cell[r].append((idx, l, allowed, False))

    def prune(cell, other, allowed, is_left):
        m = allowed if is_left else [[allowed[j][i] for j in range(3)]
                                     for i in range(3)]
        new = 0
        for d in range(3):
            if dom[cell] >> d & 1:
                if any(m[d][e] and (dom[other] >> e & 1) for e in range(3)):
                    new |= 1 << d
        return new

    changed = True
    rounds = 0
    empty_at = None
    while changed and rounds < 200:
        changed = False
        rounds += 1
        for c in range(z.n_cells):
            cur = dom[c]
            if cur == 0:
                empty_at = c
                break
            new = 0
            for d in range(3):
                if not (cur >> d & 1):
                    continue
                keep = True
                for (idx, o, allowed, is_left) in by_cell[c]:
                    if prune(c, o, allowed, is_left) & (1 << d) == 0:
                        keep = False
                        break
                if keep:
                    new |= 1 << d
            if new != cur:
                dom[c] = new
                changed = True
            if new == 0:
                empty_at = c
                break
        if empty_at is not None:
            break

    free = sum(1 for d in dom if d == 0)
    fixed = sum(1 for d in dom if bin(d).count("1") == 1)
    print("propagation rounds=%d domains_empty=%d fully_fixed=%d"
          % (rounds, free, fixed))
    if empty_at is not None:
        c = empty_at
        print("NO_ASSIGNMENT cell=%d" % c)
        for (idx, o, allowed, is_left) in by_cell[c][:4]:
            l, r, _ = cons[idx]
            print("  blocking face idx=%d pair=(%d,%d) allowed=%s"
                  % (idx, l, r, allowed))
        return 1
    # greedy completion then full verification
    import random
    random.seed(0)
    val = [0] * z.n_cells
    for c in range(z.n_cells):
        opts = [d for d in range(3) if dom[c] >> d & 1]
        val[c] = random.choice(opts)
    for _ in range(200):
        bad = 0
        for idx, (l, r, allowed) in enumerate(cons):
            if not allowed[val[l]][val[r]]:
                val[r] = random.choice([d for d in range(3)
                                        if dom[r] >> d & 1 and allowed[val[l]][d]])
                bad += 1
        if bad == 0:
            break
    viol = sum(1 for (l, r, allowed) in cons if not allowed[val[l]][val[r]])
    print("ASSIGNMENT violations=%d" % viol)
    if viol:
        print("NO_ASSIGNMENT (greedy failed; constraints may need full CSP)")
        return 1
    print("MERGE_DIRS_PASS mode=%s" % mode)
    return 0


if __name__ == "__main__":
    sys.exit(main())
