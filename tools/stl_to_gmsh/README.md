# Planar STL to Gmsh 2D geometry

Reads ASCII or binary STL in the XY plane, welds coincident vertices, extracts
one closed boundary loop, and writes a Gmsh geometry file. Gmsh is an optional
external executable; OneFLOW does not link the Gmsh SDK.

STL has no standard CFD boundary-condition metadata. The JSON sidecar assigns
every extracted edge to a named physical curve. The square fixture labels its
sides Inlet, Outlet, Symmetry and Farfield and requests one quadrilateral cell.

Run from the repository root:

    python tools/stl_to_gmsh/stl_to_gmsh.py tools/stl_to_gmsh/data/square.stl --boundaries tools/stl_to_gmsh/data/square.boundary.json --geo square.geo --gmsh gmsh --msh square.msh

Tests do not require Gmsh:

    python -m unittest tools.stl_to_gmsh.test_stl_to_gmsh

Limits: one simple closed loop in a constant-Z plane; no holes, multiple
components, arbitrary plane projection, curved-boundary reconstruction or
general multi-cell quad layout. The one-quad fixture is for topology and
conversion validation, not a useful CFD discretization by itself.
