# 2D Channel Internal Flow

This reference Case treats the closed rectangular STL loop as the boundary of
the fluid region. Its four sides are assigned to `Inlet`, `Outlet`,
`LowerWall` and `UpperWall`. A transfinite 20 x 5 layout creates a
reproducible 100-cell quadrilateral channel mesh.

Run from the repository root:

    python grid/stl_to_gmsh/run_python_case.py --case-directory grid/stl_to_gmsh/channel_internal_flow/python

The expected node, boundary-line and quadrilateral counts are exact. This is
a geometry/mesh integration Case, not yet a CFD physical-validation Case.
