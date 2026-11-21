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
* The plane points and normals are explicit inputs in this slice. Automatic
  local direction estimation and interactive 3-D endpoint selection are
  subsequent GUI work.
* This stage produces geometry only. Gmsh volume meshing and the OneFLOW CFD
  smoke run remain manual Mesh Generation / Mesh-to-CFD Validation actions;
  they are not CI jobs.

Verification
------------

Run the pure geometry regression tests without Gmsh:

.. code-block:: console

   python -m unittest grid.tests.python.stl_import.test_generic_stl_surface_cut

The next integration slice will write the cut-and-capped STL to the user's
.oneflow workspace, preserve the three boundary identities, and hand the
closed surface to the existing Gmsh volume-mesh workflow.


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

The GUI can open the closed surface in the existing STL-to-Gmsh page. This
hands off a closed surface for manual volume meshing; preserving inlet, outlet
and wall as distinct Gmsh Physical Surface groups is a subsequent integration
step. The current Gmsh volume-mesh check remains manual and is not part of CI.
