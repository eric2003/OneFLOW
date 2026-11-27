# Square-duct T-junction

A controlled planar counterpart to the circular OECD/NEA–Vattenfall T-junction.
The main duct is 4 m long with a 1 m × 1 m cross-section; the branch is 2 m
long with the same cross-section and joins at 90 degrees.

The STL is built from explicit coordinates and triangle connectivity. No OCC
Boolean intersection or curved-surface intersection is used. Shared vertices
are emitted from identical coordinates. The generator checks closed-shell edge
incidence and positive signed volume, and writes checksum-bound region and
three-port subregion manifests.

Operating modes: `main_inlets` (both main ends inflow, branch outflow),
`branch_inlet` (branch inflow, both main ends outflow), and `benchmark`
(left main and branch inflow, right main outflow). The initial Gmsh target is
deliberately coarse (0.75 m) and uses TETRA_4. Windows CI runs a two-step CFD
smoke test; this is not physical validation.
