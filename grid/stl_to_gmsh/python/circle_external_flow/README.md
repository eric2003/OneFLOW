# 2D Circle External Flow

This reference Case uses a planar, 32-segment circular STL obstacle inside a
rectangular outer domain. The obstacle loop becomes the internal `Wall`
boundary; the outer rectangle is split into `Inlet`, `Outlet`,
`Symmetry` and `Farfield` Physical Curves. The remaining surface is the
`Fluid` region.

Run from the repository root:

    python grid/stl_to_gmsh/python/run_case.py --case-directory grid/stl_to_gmsh/python/circle_external_flow

The generated unstructured mesh is validated by Physical Group names and
minimum node/line/triangle counts. Exact counts are intentionally not pinned
because Gmsh's unstructured meshing details may vary by release.
