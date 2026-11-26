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
text. A separate opt-in integration test invokes Gmsh 4.15.2 on both the blind
cylinder cut and a non-cylindrical prolate shell cut at two planes. The latter
has a curved wall with two boundary loops. Both tests check TETRA_4 volume
cells and exterior TRI_3 facets, verify the four Physical Groups, and require
every exterior triangle to belong to exactly one of `wall`, `inflow` or
`outflow`. The test runs on both Linux and Windows
in GitHub Actions. CGNS inspection and OneFLOW CFD boundary verification
remain downstream checks; a successful solver launch alone is not accepted
as evidence that the three physical boundary conditions are correct.


Oblique cylinder inlet reference case
-------------------------------------

The standard reference case is stored under
``grid/stl_to_gmsh/3d/cylinder_internal_flow/oblique_cut/python``.
It reuses the canonical closed cylinder STL as a blind triangle-shell input;
there is deliberately no geometry sidecar for the cutting operation.

The cylinder axis is X and its reference radius is ``R = 1``. The inlet
plane passes through ``(-1.2, 0, 0)`` with unit normal
``(cos(30 deg), sin(30 deg), 0)``. The outlet remains at ``x = 1.5``
with normal ``(1, 0, 0)``.

For the inlet plane,

``cos(theta) (x - x0) + sin(theta) y = 0``

and the cylinder satisfies ``y^2 + z^2 = R^2``. In plane-local
coordinates, the ellipse semi-axes are ``R / cos(theta)`` and ``R``.
At 30 degrees the expected major/minor ratio is approximately 1.1547. The
outlet remains circular. These analytic values are used only by verification;
the cutter receives only STL triangles and plane points/normals.

Run the complete cut and Gmsh volume-mesh example from the repository root:

.. code-block:: console

   python grid/stl_to_gmsh/3d/cylinder_internal_flow/oblique_cut/python/run_cylinder_oblique_cut.py --gmsh gmsh

The runner writes generated geometry and mesh artifacts below the user's
``.oneflow/mesh-generation/cylinder_oblique_cut`` directory. It retains
``wall``, ``inflow`` and ``outflow`` as separate named surface regions
and requires the ``fluid_domain`` volume group. The pure geometry test checks
the measured inlet ellipse and circular outlet; the opt-in Gmsh integration
checks non-empty Physical Groups and complete exterior-facet assignment.

The Mesh Generation GUI now offers **Convert + Run OneFLOW** after this MSH is
generated. It hands the existing mesh to the shared Mesh-to-CFD Validation
pipeline, which stages the mesh, converts Gmsh MSH to CGNS, converts CGNS to
the native OneFLOW grid, and runs a short 3-D solver smoke test. Validation
artifacts and the run report are written below the mesh project's ``validation``
directory. The smoke test confirms workflow execution and mesh ingestion only;
it is not a claim of physically validated oblique-inlet CFD accuracy.


Native C++ surface-cut entry
----------------------------

The native geometry operation is implemented under ``codes/geometry/stl``. STL cutting is opt-in through the grid control parameters: ``stlSurfaceCutEnabled = 1`` enables the operation, ``stlSurfaceCutConfigFile`` identifies a JSON file relative to the case, and ``stlSurfaceCutTolerance`` supplies the fallback geometric tolerance. The JSON uses ``cutting_planes.inlet/outlet.point`` and ``normal``, with ``keep_positive`` selecting the retained half-space. The default is disabled, so existing STL workflows are unchanged. The current native implementation requires one connected, closed, outward-oriented input shell and one simple loop at each cut plane.


Native Case Execution
---------------------

The native case is registered as
``stl_to_gmsh_cylinder_internal_flow_oblique_cut`` in
``grid/case_catalog.json``. It is a normal OneFLOW C++ Case under
``grid/stl_to_gmsh/3d/cylinder_internal_flow/oblique_cut/cpp``, not a
standalone utility. It can be selected in Run Development Cases or directly
in Mesh-to-CFD Validation. The validation workflow stages the Case in the
user's OneFLOW workspace, stages the shared parent-level cut configuration,
runs the native mesh producer, checks the four Physical Groups, then continues
through MSH → CGNS → OneFLOW grid → CFD smoke test.

The Python reference entry is separately registered as
``stl_to_gmsh_python_cylinder_internal_flow_oblique_cut``. Its runner is
``grid/stl_to_gmsh/3d/cylinder_internal_flow/oblique_cut/python/run_cylinder_oblique_cut.py``.
Both implementations read the common
``grid/stl_to_gmsh/3d/cylinder_internal_flow/oblique_cut/script/cutting.json``.

The native OneFLOW parameters are:

* ``stlSurfaceCutEnabled`` (default ``0``): set to 1 to enable native
  plane cutting.
* ``stlSurfaceCutConfigFile`` (default empty): Case-relative path to the
  shared cutting JSON.
* ``stlSurfaceCutTolerance`` (default ``1e-8``): fallback tolerance;
  an explicit JSON tolerance takes precedence.

When disabled, the existing STL conversion behavior is unchanged. When
enabled, the current native cutter requires a single connected, closed,
outward-oriented shell and one simple boundary loop at each cutting plane.
The native operation requires a 3-D layout with
``volume_mesh = closed_surface_automatic_tetrahedra``. The current Gmsh volume
converter additionally requires genus-zero input and rejects zero-area or
near-degenerate triangles.
