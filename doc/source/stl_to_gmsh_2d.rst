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

grid/stl currently provides:

* ASCII and binary STL readers;
* scale-controlled vertex welding and planar XY validation;
* triangle-edge incidence analysis and closed-loop boundary reconstruction;
* rejection of open/branched boundaries, multiple loops, non-manifold edges,
  zero-area triangles and self-intersecting outer polygons;
* a sidecar reader that requires every outer edge to be assigned exactly once;
* Gmsh geometry generation, including physical curves and an optional single
  quadrilateral cell configuration;
* optional invocation of a separately installed Gmsh executable.

External-flow obstacle domain
------------------------------

Use a sidecar with explicit outer bounds and one name for each rectangle side:

.. code-block:: json

   {
     "mesh_size": 0.25,
     "surface_name": "Fluid",
     "body_boundary_name": "Wall",
     "domain": {
       "xmin": -1.0, "xmax": 2.0,
       "ymin": -1.0, "ymax": 2.0,
       "boundaries": {
         "left": "Inlet",
         "right": "Outlet",
         "bottom": "Symmetry",
         "top": "Farfield"
       }
     }
   }

The rectangle must strictly enclose the STL body. Its four sides must have
distinct, non-empty physical names. The generated surface has the STL loop as
a hole, so the body is a wall boundary and the outer rectangle bounds the
fluid region. This mode generates an unstructured 2D mesh; the one-quad option
remains available only for the original four-edge square fixture.

Example
-------

.. code-block:: bash

   python grid/stl/stl_to_gmsh.py grid/stl/fixtures/square.stl --boundaries grid/stl/fixtures/square.boundary.json --geo /tmp/square.geo --gmsh gmsh --msh /tmp/square.msh

Parser and geometry-generation tests do not require Gmsh:

.. code-block:: bash

   python -m unittest tools.stl_to_gmsh.test_stl_to_gmsh

Current limits and next steps
-----------------------------

The implementation accepts one simple closed STL loop in a constant-Z plane.
The external-flow mode constructs a rectangular outer boundary with explicit
side names and treats the STL loop as one internal wall. Arbitrarily oriented
planar projection, multiple body loops, holes inside the STL body, curved
boundary reconstruction, automatic farfield sizing, mesh-size fields and
general multi-cell structured layouts are not yet implemented. The one-quad
option is a topology and conversion fixture, not a useful CFD discretization
by itself.

The next integration step is to run the generated MSH through the existing
Gmsh MSH reader and shared unstructured CGNS writer, then validate that the
four physical curve groups survive as four CGNS boundary patches and are
mapped into the intended OneFLOW boundary records.
