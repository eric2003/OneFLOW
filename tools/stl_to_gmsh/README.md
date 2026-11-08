# Planar STL to Gmsh 2D geometry

Reads ASCII or binary STL in the XY plane, welds coincident vertices, extracts
one closed boundary loop, and writes a Gmsh geometry file. Gmsh is an optional
external executable; OneFLOW does not link the Gmsh SDK.

STL has no standard CFD boundary-condition metadata. The JSON sidecar either
assigns each edge of a planar region to a physical curve, or defines an
external-flow rectangle around an STL body. In external-flow mode the body is
an internal `Wall` hole and the rectangle sides are named `Inlet`, `Outlet`,
`Symmetry` and `Farfield`.

Run from the repository root:

    python tools/stl_to_gmsh/stl_to_gmsh.py tools/stl_to_gmsh/data/square.stl --boundaries tools/stl_to_gmsh/data/square.boundary.json --geo square.geo --gmsh gmsh --msh square.msh

Tests do not require Gmsh:

    python -m unittest tools.stl_to_gmsh.test_stl_to_gmsh

External-flow sidecar example:

    {
      "mesh_size": 0.25,
      "domain": {
        "xmin": -1, "xmax": 2, "ymin": -1, "ymax": 2,
        "boundaries": {
          "left": "Inlet", "right": "Outlet",
          "bottom": "Symmetry", "top": "Farfield"
        }
      },
      "body_boundary_name": "Wall"
    }

The rectangle must strictly enclose the STL body. Limits: one simple closed
loop in a constant-Z plane; no multiple body loops, arbitrary plane projection,
curved-boundary reconstruction or general multi-cell structured layout. The
one-quad fixture is for topology and conversion validation, not a useful CFD
discretization by itself.
