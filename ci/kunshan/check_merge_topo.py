"""Topology feasibility checker for a per-parent local merge of a uniformly
refined (2x2x2) .ofl mesh.  Writes no grid.

Provenance rule under test (refine_ofl.py): fine cell id = parent*8 + local
corner k, where k indexes the parent frame corner in CHILD_LATTICE order.
If that rule holds and neighbouring parents share the same axis convention,
then any two fine cells that share an interior face must have octant
coordinates differing in exactly one component.

  same parent (l//8 == r//8) : must differ in exactly 1  (intra-cell cut)
  different parents          : must differ in exactly 1  (frames aligned)
                               otherwise the two parents disagree about which
                               axis crosses their shared face -> a per-parent
                               local merge pattern cannot be compatible there.

Every offending face is printed as a minimal counterexample.
"""
import sys
from collections import Counter
sys.path.insert(0, "ci/kunshan")
import refine_ofl as R
import agglomerate_ofl as A

CH = A.CHILD_LATTICE


def ndiff(a, b):
    return sum(1 for c in range(3) if CH[a][c] != CH[b][c])


def main(path):
    z, _, _ = R.parse_ofl(path)
    if z.n_cells % 8:
        print("FAIL subtle: cells=%d not a multiple of 8" % z.n_cells)
        return 1
    intra = Counter()
    cross = Counter()
    bad_intra, bad_cross = [], []
    for f in range(z.n_bfaces, z.n_faces):
        l, r = z.left[f], z.right[f]
        n = ndiff(l % 8, r % 8)
        if l // 8 == r // 8:
            intra[n] += 1
            if n != 1 and len(bad_intra) < 5:
                bad_intra.append((f, l, r, l % 8, r % 8, n))
        else:
            cross[n] += 1
            if n != 1 and len(bad_cross) < 5:
                bad_cross.append((f, l, r, l // 8, r // 8, l % 8, r % 8, n))
    print("cells=%d parents=%d faces=%d nBFaces=%d"
          % (z.n_cells, z.n_cells // 8, z.n_faces, z.n_bfaces))
    print("intra-parent face octant-difference histogram: %s" % dict(intra))
    print("cross-parent face octant-difference histogram: %s" % dict(cross))
    for b in bad_intra:
        print("COUNTEREXAMPLE intra face=%d fine=(%d,%d) octants=(%d,%d) diff=%d" % b)
    for b in bad_cross:
        print("COUNTEREXAMPLE cross face=%d fine=(%d,%d) parents=(%d,%d) "
              "octants=(%d,%d) diff=%d" % b)
    ok = not bad_intra and not bad_cross
    print("MERGE_TOPO_%s" % ("PASS" if ok else "CONFLICT"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1]))
