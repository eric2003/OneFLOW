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

Native OneFLOW Cases
--------------------

The primary runtime implementation is now available in two normal OneFLOW
Grid Cases:

* `grid/stl_to_gmsh/square/single_quad/cpp`: one square, one QUAD_4 cell, four boundary edges;
* `grid/stl_to_gmsh/square/two_triangles/cpp`: the same square, two TRI_3 cells;
* `grid/stl_to_gmsh/two_square_cells/cpp`: a 2 x 1 rectangle, two QUAD_4 cells,
  six exterior boundary edges.

Both Cases use the standard C++ STL reader and Grid conversion path. The
Case's `script/grid.txt` selects `sourceGridType = "stl"`,
`targetGridType = "gmsh"`, and `gmshExecutable`. The JSON
`gridLayoutFileName` sidecar defines `transfinite_cells` and assigns the
four rectangle sides to `inflow`, `outflow`, `symmetry`, and `farfield`.
The fluid surface is named `fluid_domain`.

Run native and Python Cases only through the external staging orchestrator
after building OneFLOW and installing Gmsh:

.. code-block:: console

   python grid/stl_to_gmsh/develop_test/run_cases.py --oneflow-executable build/bin/Debug/OneFLOW.exe --gmsh gmsh

The orchestrator copies each Case to the user's .oneflow staging directory
and launches OneFLOW with that staged Case as its working directory. Do not
run checked-in Case directories directly: generated files must remain outside
the repository source tree.

The native path writes .geo, MSH 4.1 and a validation report under each
staged Case's results directory. It reads the generated MSH again with the existing
C++ Gmsh reader and checks node, boundary-line, QUAD_4 and Physical Group
counts. Independent Python reference Cases are maintained under
`grid/stl_to_gmsh/<case>/python/`; the existing external-flow Python Case remains at
`grid/stl_to_gmsh/script/run.py` with `config/external_flow.json`.

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
a Physical Surface group named ``fluid_domain``. Physical names are mesh metadata, not
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
* Gmsh geometry generation, including physical curves, transfinite structured
  layouts, and automatic unstructured triangle meshing;
* optional invocation of a separately installed Gmsh executable.

Automatic triangular meshing
-----------------------------

For a directly bounded planar region, omit the transfinite_cells field from the
layout. The converter then leaves the surface unstructured and lets Gmsh
choose the triangulation. mesh_size controls the target characteristic length;
it is not a requested element count. Smaller values normally produce more,
smaller elements, but exact node and triangle counts depend on the Gmsh version,
meshing algorithm and geometry.

The automatic mode sets Mesh.MeshSizeMin and Mesh.MeshSizeMax to the configured
value, enables point-based sizing and boundary-size extension, and does not
apply transfinite constraints or recombination. C++ and Python use the same
layout semantics. Validation checks non-empty first-order triangles, boundary
LINE_2 correspondence, manifold edge incidence, the disk Euler characteristic,
and declared Physical Groups rather than exact element counts.

Use transfinite_cells only when a specific structured topology is required.
In that mode integer x/y cell counts control the structured layout and
recombine_surface selects QUAD_4 versus TRI_3 output.

External-flow obstacle domain
------------------------------

Use a sidecar with explicit outer bounds and one name for each rectangle side:

.. code-block:: json

   {
     "mesh_size": 0.25,
     "surface_name": "fluid_domain",
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

   python -m unittest grid.stl.test_stl_to_gmsh

Native C++ circle external-flow Case
------------------------------------

The native Case at `grid/stl_to_gmsh/circle_external_flow/cpp` uses the same
32-segment circular STL and layout as the independent Python reference at
`grid/stl_to_gmsh/circle_external_flow/python`. Its layout defines an outer
rectangular computational domain and the four outer boundary labels. The C++
converter creates the outer loop, reverses the body loop orientation to form
a hole, and creates one fluid Plane Surface. The upper and lower outer curves
share the `farfield` Physical Group; the body loop is `wall`.

The native Case is registered as a C++ development Case and is automatically
available in Mesh-to-CFD Validation. The validation graph uses the selected
Case's declared MSH path for the Gmsh-to-CGNS and Gmsh-to-OneFLOW handoffs.
The default CFD consumer is a solver smoke test only; a physical cylinder
wake validation requires a separate case setup and reference quantities.

The C++ and Python Case directories intentionally mirror one another. Their
`input/model.stl` files are byte-identical and their `script/layout.json`
files are identical. The layout is the shared geometry contract: it defines
the body boundary name, outer rectangle, side-to-Physical-Group assignments,
surface name, body/outer mesh sizes, and quality limits. The mesh comparison
is semantic rather than byte-for-byte: Gmsh versions and independently
generated .geo files may produce different valid triangulations.

Both implementations verify that generated MSH 4.1 coordinates reproduce
the configured domain bounds and remain planar in XY. They require the
`wall`, `inflow`, `outflow`, `farfield`, and `fluid_domain` Physical
Groups with the expected dimensions, and reject degenerate TRI_3 cells.
Triangle quality uses the normalized measure

.. math::

   q = \frac{4\sqrt{3} A}{a^2+b^2+c^2}

where :math:`A` is triangle area and :math:`a,b,c` are its edge lengths.
An equilateral triangle has :math:`q=1`; values approach zero for
degenerate triangles. The shared case contract requires :math:`q \ge 0.01`
and a minimum interior angle of at least one degree. The Python MSH inspector
reports minimum area, minimum normalized quality, minimum angle, coordinate
bounds and orientation counts. The native C++ path applies the same quality
limits and also checks its mesh-read report for degenerate cells.

External- and internal-flow reference Cases
--------------------------------------------

Two independent Python reference Cases make the computational-domain
distinction executable:

* `grid/stl_to_gmsh/circle_external_flow/cpp` and `grid/stl_to_gmsh/circle_external_flow/python`:
  a 32-segment circular STL obstacle is placed inside a rectangular outer
  domain. The fluid surface is the outer loop with the obstacle loop as a
  hole. The body is `wall`; the outer curves are `inflow`, `outflow` and
  `farfield` (the top and bottom edges share the farfield group).
* `grid/stl_to_gmsh/channel_internal_flow/python`: a rectangular STL loop
  directly bounds the fluid region. Its four edges are `Inlet`, `Outlet`,
  `LowerWall` and `UpperWall`. A transfinite 20 x 5 layout produces 100
  quadrilateral cells.

Both Cases use the same STL reader, boundary reconstruction, sidecar
validation, Gmsh geometry writer and MSH 4.1 inspector. The channel Case has
exact expected mesh counts. The external unstructured Case validates physical
groups and minimum non-zero mesh counts instead of pinning an exact count that
may vary with Gmsh releases.

These are geometry and mesh-generation reference Cases, not complete CFD
validation cases. They do not yet establish boundary-condition state mapping,
solver convergence, or physical accuracy. The next integration step is to
pass their generated MSH artifacts through the existing Gmsh reader and CGNS
writer, then validate that physical curve groups survive as distinct CGNS
boundary patches and are mapped into intended OneFLOW boundary records.

Current limits and next steps
-----------------------------

The 2-D implementation accepts one simple closed STL loop in a constant-Z
plane. The external-flow mode constructs a rectangular outer boundary with
explicit side names and treats the STL loop as one internal wall. Arbitrarily
oriented planar projection, multiple body loops, holes inside the STL body,
curved boundary reconstruction, automatic farfield sizing, mesh-size fields
and general multi-cell structured layouts are not yet implemented.

A separate 3-D `cube_single_hex` prototype now accepts only a closed,
axis-aligned cube STL with equal edge lengths, exactly 12 triangles, eight
bounding-box corners, two triangles on each face, and closed two-manifold
edge incidence. It reconstructs six planar surfaces, creates one volume, and
uses Gmsh transfinite/recombined meshing. The exact acceptance topology is
eight nodes, six QUAD_4 boundary faces and one HEXA_8 volume cell. XMin,
XMax, YMin and YMax map to inflow, outflow, symmetry and farfield; ZMin and
ZMax share wall; fluid_domain names the volume.

The cube is a restricted geometry-recovery prototype, not a general STL
volume mesher. Mesh-to-CFD Validation passes it through MSH, CGNS and OneFLOW
grid conversion and verifies FaceCenter boundary semantics and the boundary
region map. The current validation graph then runs a short, dimension-adjusted
solver smoke test using shared consumer parameters. A single-cell cube remains
a topology/interoperability fixture, not a meaningful CFD accuracy case.


Automatic square triangulation Case
-----------------------------------

The square automatic-triangulation native C++ Case and independent Python
reference Case are registered in the shared Case catalog. Both reuse the
same two-facet square STL and four outer Physical Curve groups. They omit
``transfinite_cells`` and ``recombine_surface``. Gmsh determines the TRI_3
topology from the configured ``mesh_size``, which is a target characteristic
length rather than an exact triangle count. Both implementations are
selectable through Run Development Cases / Python reference cases and either
can be selected in Mesh-to-CFD Validation. The downstream CFD consumer remains
a smoke test, not an accuracy benchmark for this geometry.
