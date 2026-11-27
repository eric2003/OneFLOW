# Square-duct T-junction

A controlled planar counterpart to the circular OECD/NEA–Vattenfall T-junction.
The main duct is 4 m long with a 1 m × 1 m cross-section; the branch is 2 m
long with the same cross-section and joins at 90 degrees.

The STL is built from explicit coordinates and triangle connectivity. No OCC
Boolean intersection or curved-surface intersection is used. Shared vertices
are emitted from identical coordinates. The generator checks closed-shell edge
incidence and positive signed volume, and writes checksum-bound region and
three-port subregion manifests.

## Operating modes

- `main_inlets`: both main-duct ends are inflows; branch is outflow.
- `branch_inlet`: branch is inflow; both main-duct ends are outflows.
- `benchmark`: left main port and branch are inflows; right main port is outflow.

## Run

```powershell
python grid/stl_to_gmsh/3d/pipe_tee/square_duct/automatic_tetrahedra/run_square_tee_case.py --case-directory "$env:USERPROFILE/.oneflow/mesh-validation/square-duct-tee" --gmsh gmsh --mesh-size 0.75 --operating-mode main_inlets
```

The initial target is deliberately coarse and produces a TETRA_4 mesh. CFD
scope is a bounded two-step smoke test only, not physical validation. Generated
artifacts belong in the OneFLOW user-data directory, not this source tree.
