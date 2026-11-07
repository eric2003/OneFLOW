Planar STL to Gmsh 2D Geometry
==============================

Purpose
-------

The first STL-to-Gmsh milestone is deliberately two-dimensional. It reads a
triangulated STL region in the XY plane, reconstructs its single closed outer
boundary, and writes a Gmsh .geo geometry. Gmsh remains an external optional
executable; OneFLOW does not link the Gmsh SDK for this workflow.

The first fixture is a unit square represented by two STL triangles. Their
shared diagonal is an interior edge; the four singly incident edges form the
outer boundary. The generated Gmsh geometry can request one quadrilateral CELL
using transfinite curves, a transfinite surface, and recombination.

Input and boundary semantics
----------------------------

STL stores triangular facets and does not define a standard CFD boundary
condition vocabulary. It cannot, by itself, reliably say that one side is an
inlet and another is an outlet. The workflow therefore uses two inputs:

* square.stl: coordinates and triangulated surface topology;
* square.boundary.json: explicit named assignments for every extracted
  boundary edge, plus meshing options.

The example assigns the four square edges to Inlet, Outlet, Symmetry and
Farfield. These become Gmsh Physical Curve groups. The 2D fluid region becomes
a Physical Surface group named Fluid. Physical names are mesh metadata, not
solver boundary-condition types. A subsequent explicit mapping layer must map
these names to OneFLOW boundary semantics and CGNS BC types.

Implementation
--------------

tools/stl_to_gmsh currently provides:

* ASCII and binary STL readers;
* scale-controlled vertex welding and planar XY validation;
* triangle-edge incidence analysis and closed-loop boundary reconstruction;
* rejection of open/branched boundaries, multiple loops, non-manifold edges,
  zero-area triangles and self-intersecting outer polygons;
* a sidecar reader that requires every outer edge to be assigned exactly once;
* Gmsh geometry generation, including physical curves and an optional single
  quadrilateral cell configuration;
* optional invocation of a separately installed Gmsh executable.

Example
-------

.. code-block:: bash

   python tools/stl_to_gmsh/stl_to_gmsh.py tools/stl_to_gmsh/data/square.stl --boundaries tools/stl_to_gmsh/data/square.boundary.json --geo /tmp/square.geo --gmsh gmsh --msh /tmp/square.msh

Parser and geometry-generation tests do not require Gmsh:

.. code-block:: bash

   python -m unittest tools.stl_to_gmsh.test_stl_to_gmsh

Current limits and next steps
-----------------------------

The initial implementation accepts one simple closed loop in a constant-Z
plane. Arbitrarily oriented planar projection, holes, multiple disconnected
regions, curved boundary reconstruction, automatic farfield construction,
mesh-size fields and general multi-cell quadrilateral layouts are not yet
implemented. The one-quad option is a topology and conversion fixture, not a
useful CFD discretization by itself.

The next integration step is to run the generated MSH through the existing
Gmsh MSH reader and shared unstructured CGNS writer, then validate that the
four physical curve groups survive as four CGNS boundary patches and are
mapped into the intended OneFLOW boundary records.
