Generic STL Surface Cutting
=============================

Purpose
-------

This stage starts with the canonical internal-flow cylinder STL, but treats it
as an unknown triangulated surface. The cutter must not read the neighboring
.geometry.json manifest, use cylinder parameters, or evaluate an analytic
cylinder equation. The cylinder is a convenient first fixture, not a special
geometry type.

The implementation is in grid/stl_import/surface_cut.py. It is a geometry
operation over STL triangles and two user-supplied planes.

Input contract
--------------

* One connected, closed, consistently oriented, two-manifold triangle shell.
* An inlet point and a unit-normal direction. The inlet normal points into the
  retained fluid domain.
* An outlet point and a unit-normal direction. The outlet normal points from
  the retained domain toward the outlet.
* A positive geometric tolerance used for vertex welding and plane tests.

Cutting algorithm
-----------------

For each source triangle, the implementation clips its polygon against the
positive half-space of the inlet plane and the negative half-space of the
outlet plane. Intersections are computed from signed point-to-plane distances;
no primitive-specific formula is used. Resulting polygons are triangulated
without changing their source winding.

The one-use edges of the clipped side surface are collected and classified by
their distance to the two cutting planes. Each plane must produce exactly one
closed, non-branching boundary loop in this first implementation. The loop is
projected into a local two-dimensional basis and triangulated with ear clipping.
The inlet cap is oriented outward opposite the inlet normal; the outlet cap is
oriented along the outlet normal.

The returned data includes the clipped side triangles, ordered inlet/outlet
loops, cap triangles and the combined closed surface. The caller can pass the
closed surface through the existing closed-shell topology validation before
Gmsh volume meshing.

Current scope and limitations
-----------------------------

* The first reference input is the canonical cylinder STL under
  grid/stl_to_gmsh/3d/cylinder_internal_flow/automatic_tetrahedra/python/input.
* The test deliberately reads only model.stl. It does not read the analytic
  geometry manifest or patch-triangle indices.
* One simple loop per cutting plane is supported. Branched pipes, multiple
  loops, holes in a cross-section and self-intersecting cuts are rejected.
* The plane points and normals are explicit GUI inputs. Automatic local
  direction estimation and picking endpoints directly on a 3-D model are not
  included.
* The cutter produces geometry only. The generic STL-to-Gmsh volume-mesh
  path is exercised by an opt-in integration test on Linux and Windows using
  Gmsh 4.15.2. CGNS conversion and the OneFLOW CFD smoke run remain downstream
  Mesh-to-CFD Validation stages.
* The separate wall/inlet/outlet STL artifacts retain the three surface
  regions. The project importer carries their checksummed manifest and patch
  files into the mesh project, where the region-aware converter assigns
  distinct Gmsh Physical Surface groups.

Verification
------------

Run the pure geometry regression tests without Gmsh:

.. code-block:: console

   python -m unittest grid.tests.python.stl_import.test_generic_stl_surface_cut

Run the same suite with the Gmsh integration test enabled (Gmsh 4.15.2 must
be available on ``PATH``):

.. code-block:: console

   ONEFLOW_RUN_GMSH_TESTS=1 python -m unittest grid.tests.python.stl_import.test_generic_stl_surface_cut

The current CLI and GUI write the cut-and-capped STL and the three separate
surface-region files to the user's .oneflow workspace. The closed surface can
be opened in the existing Gmsh volume-mesh workflow. Distinct Physical Surface
groups and automatic CGNS boundary-name transfer are the next integration
steps.


Workbench operation
-------------------

The Mesh Generation hub exposes **Generic STL Surface Cutting** separately
from the analytic-cylinder demonstration. Select any input STL and enter the
inlet/outlet plane points and normals. The command preview shows the exact
Python invocation. Generated artifacts are written under the user's
`.oneflow/mesh-generation/generic_surface_cut` workspace by default:

* `wall.stl`: clipped lateral surface;
* `inlet.stl` and `outlet.stl`: generated end caps;
* `closed_surface.stl`: combined closed shell for the existing STL-to-Gmsh
  workflow;
* `cut_report.json`: inputs, counts, output paths and explicit confirmation
  that analytic geometry metadata was not read.

The GUI can open the closed surface in the existing STL-to-Gmsh page. The
project importer carries the checksum-verified region manifest and patch STL
files into the mesh project. The region-aware converter first merges coincident patch nodes, reconstructs
shared discrete topology with `CreateTopology`, and creates the geometric
surfaces with `CreateGeometry`. It then creates a closed surface loop and
volume with distinct `wall`, `inflow` and `outflow` Physical Surface groups
plus the `fluid_domain` Physical Volume. The end-to-end Gmsh path is
covered by the dedicated Linux and Windows integration job.


Boundary identity handoff
-------------------------

The cut CLI writes `region_manifest.json` alongside the four STL artifacts.
The manifest uses the schema `oneflow.stl-surface-regions` and records the
wall, inflow and outflow STL file names, triangle counts and SHA-256 hashes,
together with the closed-shell hash. The analytic geometry manifest is not
consulted.

When a user opens `closed_surface.stl` in STL to Gmsh, the mesh-project
importer copies the region manifest and its three patch STL files beside the
project's `input/model.stl`. The generic region-aware meshing path imports
the patches separately and assigns the Gmsh Physical Surface names
`wall`, `inflow` and `outflow`; the volume is `fluid_domain`. The existing
Gmsh-to-CGNS converter maps those names to `BCWall`, `BCInflow` and
`BCOutflow`, respectively.

The pure Python regression verifies the manifest and generated Gmsh geometry
text. A separate opt-in integration test invokes Gmsh 4.15.2, checks that
TETRA_4 volume cells and exterior TRI_3 facets exist, verifies the four
Physical Groups, and requires every exterior triangle to belong to exactly
one of `wall`, `inflow` or `outflow`. The test runs on both Linux and Windows
in GitHub Actions. CGNS inspection and OneFLOW CFD boundary verification
remain downstream checks; a successful solver launch alone is not accepted
as evidence that the three physical boundary conditions are correct.
