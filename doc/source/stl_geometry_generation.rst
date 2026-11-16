STL Geometry Generation for 3-D Mesh Cases
=============================================

Purpose
-------

This document records the geometric definitions and surface-construction
rules for the built-in curved 3-D STL geometries. It is separate from the
Gmsh volume-meshing contract.

The built-in geometry library currently provides three parameterized surface
generators:

* Sphere: a closed spherical surface.
* Prolate spheroid: a closed elongated ellipsoidal surface.
* Cylinder internal flow: a closed cylindrical fluid-domain boundary with
  separately identified inlet, outlet and wall patches.

The GUI entry is Mesh Generation -> Generate STL Geometry. Generated STL files
and their companion .geometry.json manifests are written to the user workspace,
not to the source tree. This geometry library is distinct from the checked-in
STL Case database used by STL to Gmsh; either can supply the input surface for
a mesh project. The cylinder manifest preserves patch identity, but the generic
closed-shell Gmsh adapter does not yet reconstruct those patch identities as
separate Physical Surfaces.

Geometry library defaults
-------------------------

The Generate STL Geometry dialog uses these starting values. They are
convenient reference defaults, not physical-unit requirements or fixed mesh
resolutions.

+----------------------+--------------------------+------------------------------+
| Geometry             | Analytic parameters      | Surface discretization      |
+======================+==========================+==============================+
| Sphere               | radius = 1               | icosphere subdivisions = 3   |
+----------------------+--------------------------+------------------------------+
| Prolate spheroid     | a = 3, b = 0.5, c = 0.5  | icosphere subdivisions = 3   |
+----------------------+--------------------------+------------------------------+
| Cylinder internal    | radius = 1, length = 4   | N_theta = 48, N_x = 16       |
| flow                 |                          |                              |
+----------------------+--------------------------+------------------------------+

The default coordinate unit is metres. The GUI allows m, mm, cm and in as
metadata labels; it does not silently scale coordinates when the unit changes.
The cylinder axis is the x-axis. Its radius and length are analytic dimensions,
while N_theta and N_x control only the triangulation of its boundary.

Generated assets are written by default to
~/.oneflow/mesh-generation/generated_stl/<geometry-id>/<geometry-id>.stl,
with a sibling .geometry.json manifest. ONEFLOW_HOME relocates the workspace.
The manifest records the intended Gmsh Physical Group names; it does not mean
that STL itself contains Gmsh entities or that a Gmsh volume mesh has already
been created.

Geometry parameters define the STL surface. Gmsh mesh-size parameters define
the volume mesh. Keep these controls independent so mesh-refinement studies
can reuse exactly the same STL input.

Current asset provenance
------------------------

The checked-in sphere asset is a watertight, outward-oriented triangulated
unit sphere centered at the origin with radius one and 1,280 triangular
facets. This count is consistent with an icosahedron whose faces undergo
three rounds of four-way subdivision, with new vertices projected onto the
unit sphere. The checked-in 6:1 prolate spheroid uses the same surface
topology transformed by semi-axis scale factors (3, 0.5, 0.5).

The native C++ and Python reference Cases consume byte-identical copies of
each checked-in STL. The original historical program that emitted these binary
STL files is not currently checked in as a dedicated canonical generator. The
new deterministic generator below is a reproducible replacement, not a claim
about the original fixture bytes. The icosphere construction is consistent
with the checked-in facet count and defines the new generator algorithm.

The current deterministic generator is
grid/stl_to_gmsh/generate_stl_geometry.py (version 1.0). It uses only the Python
standard library, writes binary STL, and emits a companion manifest containing
parameters, units, geometric metrics, patch-to-facet indices and the STL
SHA-256. It writes to the user workspace by default and does not overwrite
checked-in Case fixtures during normal runs. Regression tests are in
grid/tests/python/stl_to_gmsh/test_stl_geometry_generation.py.

Sphere surface construction
---------------------------

The reference sphere is defined by

.. math::

   x^2 + y^2 + z^2 = R^2

The current reference uses center (0, 0, 0) and R = 1. The implemented
deterministic triangulation is an icosphere:

1. Start with a regular icosahedron.
2. For each subdivision round, split every triangle into four triangles.
3. Normalize each new edge-midpoint vertex to radius R.
4. Reuse each midpoint vertex across adjacent triangles.
5. Orient every triangle consistently outwards.
6. Write the triangular surface as STL.

Three subdivision rounds yield 20 * 4^3 = 1,280 triangular facets,
matching the current fixture's facet count. Exact vertex ordering and STL
binary bytes are not part of the mesh Case contract; geometric and
topological invariants are.

The analytic reference measures are

.. math::

   A = 4\pi R^2, \qquad V = \frac{4}{3}\pi R^3

The STL is a piecewise-planar approximation, not an exact CAD surface.

Prolate spheroid surface construction
-------------------------------------

An axis-aligned ellipsoid centered at the origin is defined by

.. math::

   \frac{x^2}{a^2} + \frac{y^2}{b^2} + \frac{z^2}{c^2} = 1

The existing 6:1 prolate spheroid uses semi-axes a = 3, b = 0.5 and
c = 0.5. Construct it by applying the affine transform

.. math::

   (x, y, z)_{ellipsoid} = (a x_s, b y_s, c z_s)

to every vertex of the unit-sphere triangulation. Keep triangle connectivity
and winding unchanged; positive scale factors preserve orientation.

Its full length is 6, maximum diameter is 1, aspect ratio is 6:1, and its
analytic volume is

.. math::

   V = \frac{4}{3}\pi abc = \pi

The transformed triangles remain planar facets approximating the curved
ellipsoid. Reusing the same connectivity makes the sphere and spheroid a
controlled pair: topology and facet count stay fixed while shape changes.

Cylinder internal-flow surface construction
-------------------------------------------

The first cylinder generator creates a straight circular cylinder whose enclosed
volume represents the fluid region, not a solid body to be meshed. Its axis
is the x-axis, with radius R and length L:

.. math::

   -\frac{L}{2} \le x \le \frac{L}{2}, \qquad y^2 + z^2 = R^2

The closed STL shell consists of three connected patches:

* wall: cylindrical lateral surface.
* inflow: planar disk at x = -L/2.
* outflow: planar disk at x = +L/2.

The end disks close the STL shell for volume reconstruction. They are
geometric closure surfaces, not solid walls. Their future CFD boundary
identities are inflow and outflow; only the lateral patch is wall. The
enclosed region is fluid_domain.

The generator uses N_theta circumferential divisions and N_x axial intervals:

1. Create N_x + 1 circular rings at evenly spaced axial positions.
2. Create N_theta shared vertices on each ring using one angular convention.
   Represent the seam as shared topology, not duplicated coincident vertices.
3. Split each lateral quad into two triangles with a consistent diagonal
   pattern and outward radial winding.
4. Triangulate each end disk using its shared boundary ring and an interior
   center vertex (or another documented planar triangulation).
5. Reuse the exact ring vertices from the wall so there are no cracks or
   T-junctions.
6. Orient the inlet cap normal toward -x and the outlet cap normal toward +x.
7. Retain a per-triangle patch label for reconstruction as distinct Gmsh
   surface entities and Physical Groups.

The shell must be closed: every undirected mesh edge is used by exactly two
triangles with opposite directed uses. It must be one connected, genus-zero
component with positive enclosed volume. The three patch labels must cover
every triangle exactly once, and cap boundary loops must coincide exactly
with the ends of the wall patch.

STL has no portable standard Physical Group mechanism. The cylinder generator
therefore preserves patch identity in an explicit companion manifest as
per-facet index lists tied to a checksum of the STL.
The region-aware Gmsh handoff verifies the STL checksum, closed
surface topology, complete non-overlapping patch-index coverage, and the
geometric location of both end-cap patches. It then creates the corresponding
analytic OpenCASCADE cylinder and assigns separate `wall`, `inflow`, and
`outflow` Physical Surfaces, plus the `fluid_domain` Physical Volume.
This avoids classifying the inlet and outlet as solid walls. The generated
volume mesh is checked for tetrahedral cells, matching exterior facets and
the exact expected Physical Group names.

The initial geometry-only mesh Case should validate analytic bounds, radius,
length and enclosed volume; watertightness, connectedness, manifold edge
incidence and winding; exactly three non-overlapping patch labels; one
positive-volume fluid_domain; tetrahedral volume cells and matching exterior
facets; and Physical Surface groups wall, inflow and outflow plus Physical
Volume fluid_domain. No CFD solver or boundary-state configuration belongs
to this first mesh-generation stage.

Units and discretization policy
-------------------------------

STL coordinates are numbers whose physical unit must be declared by the
geometry manifest. Do not silently convert between metres and millimetres.
Keep these controls separate:

* Geometry resolution: sphere subdivision level, spheroid source-surface
  resolution, or cylinder N_theta and N_x.
* Volume-mesh resolution: Gmsh mesh_size_min and mesh_size_max.

For a mesh-size study, hold the STL bytes and checksum fixed and change only
Gmsh's target characteristic length. Geometry-resolution studies are separate
and must report the changed STL approximation.

Validation and artifact policy
-------------------------------

Every generated surface should have a machine-readable report containing the
generator name and algorithm version; geometry type, parameters and units;
vertex and triangle counts; bounding box, surface area and signed enclosed
volume where applicable; connected-component count, edge incidence and
orientation checks; patch names and per-patch triangle counts; STL checksum
and output path.

Generated STL, manifests, reports and Gmsh outputs belong in the user's
OneFLOW workspace under .oneflow/mesh-generation/, not in the source tree.
The GUI defaults to generated_stl/<geometry-id>/<geometry-id>.stl and keeps
the matching .geometry.json beside it. Checked-in STL files are immutable Case inputs during normal
runs. Fixture regeneration must be an explicit, reviewed maintenance action.

Related mesh cases
------------------

* :doc:`stl_to_gmsh_3d` — 3-D closed-shell volume meshing and acceptance limits.
* :doc:`stl_to_gmsh_sphere` — current unit-sphere Case.
* :doc:`stl_to_gmsh_prolate_spheroid` — current 6:1 spheroid Case.
* :doc:`stl_to_gmsh_commands` — direct closed-STL-to-volume-mesh command.
