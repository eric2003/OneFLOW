STL Geometry Assets and STL-to-Gmsh Pipeline
================================================

Purpose
-------

This page connects the geometry asset library to the concrete mesh-generation
adapters. It explains what an STL asset represents, how the built-in curved
surfaces are parameterized, how metadata preserves identity and boundary
patches, and how native C++ and Python paths hand the surface to Gmsh.

The current reference families are:

* Unit sphere: a closed, wall-only surface.
* 6:1 prolate spheroid: a stretched unit-sphere surface, also wall-only.
* Internal-flow cylinder: a closed shell whose lateral wall, inlet cap and
  outlet cap have distinct identities.
* 90-degree pipe bend: a swept circular section with straight inlet/outlet
  extensions; the Python generator exports a closed STL and three named patch
  STLs for the generic region-aware volume route.

An STL is a triangulated boundary representation. It is not a volume mesh,
does not contain a standard Gmsh Physical Group table, and does not encode
solver boundary-state values.

Geometry asset preparation
---------------------------

The deterministic generator is
``grid/stl_to_gmsh/generate_stl_geometry.py`` (generator version 1.0).
It uses Python standard-library code to create vertices and triangle
connectivity, orient triangles outward, calculate surface metrics, write
binary STL and write a sibling ``.geometry.json` manifest. It refuses output
inside the repository source tree. The normal destination is the user's
OneFLOW workspace, under ``.oneflow/mesh-generation/generated_stl/``.

The built-in analytic definitions and defaults are:

.. list-table::
   :header-rows: 1
   :widths: 22 36 42

   * - Asset
     - Analytic definition
     - Surface controls
   * - Unit sphere
     - `x^2 + y^2 + z^2 = R^2`, default `R=1`
     - Icosphere subdivision level, default 3
   * - Prolate spheroid
     - `x^2/a^2 + y^2/b^2 + z^2/c^2 = 1`,
       default semi-axes `(3, 0.5, 0.5)`
     - Icosphere subdivision level, default 3, then axis scaling
   * - Internal-flow cylinder
     - `-L/2 <= x <= L/2`, `y^2+z^2=R^2`,
       default `R=1, L=4`
     - Circumferential divisions 48, axial intervals 16

The sphere begins with an icosahedron. Each subdivision replaces one
triangle with four, reuses shared edge midpoints, projects new vertices onto
the radius, and preserves outward winding. Three rounds yield 1,280 facets.
The spheroid uses the same unit-sphere connectivity and scales each vertex
by `(a,b,c)`; the current reference is a prolate 6:1 body. This is a
controlled geometry pair: connectivity and facet count are held fixed while
the analytic shape changes.

The cylinder generator creates `N_x+1` axial rings with `N_theta` points
per ring, triangulates each lateral quad, and closes both ends with cap
triangles. The three patch labels are `wall`, `inflow` and `outflow`.
The inlet and outlet caps are closure surfaces for a watertight shell; they
are not solid walls.

The manifest binds the STL bytes to its asset ID, analytic definition,
parameters, units, generator name/version, surface metrics, triangle count,
patch-to-triangle index lists and SHA-256. A unit label is metadata: changing
`m` to `mm` does not scale coordinates. Geometry resolution (icosphere
subdivision or cylinder divisions) is independent of Gmsh volume mesh size.

Asset generation command examples
----------------------------------

.. code-block:: console

   python grid/stl_to_gmsh/generate_stl_geometry.py --geometry sphere --output "%USERPROFILE%/.oneflow/mesh-generation/generated_stl/sphere/sphere.stl"
   python grid/stl_to_gmsh/generate_stl_geometry.py --geometry prolate_spheroid --output "%USERPROFILE%/.oneflow/mesh-generation/generated_stl/prolate_spheroid/prolate_spheroid.stl"
   python grid/stl_to_gmsh/generate_stl_geometry.py --geometry cylinder_internal_flow --output "%USERPROFILE%/.oneflow/mesh-generation/generated_stl/cylinder_internal_flow/cylinder_internal_flow.stl"

Use the generator's ``--help`` output for all parameter overrides. Generated
STL and manifest files belong in the user workspace. Version-controlled
standard assets under ``grid/stl_to_gmsh/3d/`` are immutable Case inputs;
normal GUI runs must not regenerate them in place.

STL import and topology validation
-----------------------------------

The shared Python STL reader is ``grid/stl_import/reader.py``. It detects
ASCII or binary STL and returns triangle records without pretending that
file parsing itself validates topology. The Python surface utilities in
``grid/stl_import/surface.py`` weld vertices within a configured tolerance
and extract planar boundary loops. Closed 3-D validation is implemented in
``grid/stl_to_gmsh/converter.py``: it checks finite coordinates,
non-degenerate and non-duplicate triangles, edge incidence and orientation,
one connected component, genus-zero topology and positive signed volume.

The native C++ reader is ``StlMeshReader`` in
``codes/geometry/stl/src/StlMeshIO.cpp``, declared in
``codes/geometry/stl/include/StlMeshIO.h``. The STL-to-Gmsh dispatch and
geometry generation are in ``codes/geometry/src/GridConversion.cpp``.
The native closed-shell route validates the surface before Gmsh generation;
the cylinder route additionally checks the manifest patch partition and
analytic dimensions. Do not assume the C++ route implements every Python
manifest check: in particular, checksum verification must be documented only
where the current implementation actually performs it.

STL-to-Gmsh: method selection
----------------------------

The Case layout JSON selects the geometry pipeline using
``mesh_dimension`` and, for 3-D, ``volume_mesh``. The C++ dispatcher
``ConvertStlToGmsh`` selects the corresponding conversion routine.
Python Case execution is orchestrated by
``grid/stl_to_gmsh/run_python_case.py`` and uses
``grid/stl_to_gmsh/converter.py``.

Closed sphere and spheroid route
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For the generic closed-shell volume path:

#. Read the STL and validate a single closed, outward-oriented, genus-zero
   shell.
#. Merge the STL into the Gmsh geometry model.
#. Run ``ClassifySurfaces`` using configured surface and curve angle
   thresholds.
#. Run ``CreateGeometry`` to reconstruct parametrizable surface patches.
#. Collect the reconstructed surfaces into one ``Surface Loop``.
#. Create one enclosed ``Volume``.
#. Assign a dimension-2 ``wall`` Physical Surface and a dimension-3
   ``fluid_domain`` Physical Volume.
#. Invoke Gmsh in 3-D mode and write ASCII MSH 4.1 with first-order
   TETRA_4 cells and TRI_3 exterior facets.

The native C++ path writes a ``.geo`` script, invokes the configured Gmsh
executable, reads the generated MSH 4.1 and validates tetrahedral topology
and Physical Group names. The Python path generates equivalent Gmsh geometry
through the converter and inspects the MSH output. The sphere and spheroid
Cases exercise this generic closed-shell route.

Internal-flow cylinder route
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The cylinder uses a specialized region-aware route because classifying every
cap facet as wall would destroy the inlet/outlet meaning.

The manifest must provide the analytic radius, length, circumferential and
axial divisions, plus three non-overlapping triangle-index lists whose union
covers every STL facet exactly once. The Python converter additionally
verifies the SHA-256, expected patch counts, cap planes, lateral radius and
axial bounds. The native C++ routine
``ConvertStlParametricCylinderToGmsh`` validates the closed shell,
manifest schema, dimensions, expected patch counts, complete index
assignment, cap planes and wall radius.

After validation, the adapter constructs an analytic OpenCASCADE cylinder
rather than trying to infer three CAD faces from arbitrary STL facet
classification. It creates separate Physical Surfaces ``wall``,
``inflow`` and ``outflow``, and a Physical Volume
``fluid_domain``. Gmsh then generates the tetrahedral volume mesh. This
method preserves the CFD patch intent while using the STL and manifest as
the checked-in geometry asset contract.

Curved pipe-bend route
~~~~~~~~~~~~~~~~~~~~~~

The 90-degree pipe bend does not use an analytic OCC torus or cylinder.
Its Python geometry generator writes the closed STL, separate wall/inflow/
outflow patch STLs, a checksum-bound geometry manifest and the
``region_manifest.json`` file using the existing
``oneflow.stl-surface-regions`` schema. The generic
``generate_volume_mesh.py`` entry point detects this region manifest and
selects ``region_aware_closed_surface_automatic_tetrahedra``.

The converter verifies the closed-shell topology, the closed STL checksum,
each patch STL checksum and triangle count, and compares the orientation-
independent triangle geometry across all patch STLs against the closed STL.
This proves an exact, non-overlapping patch partition rather than relying on
matching counts alone. It then merges the three patch STLs, runs ``Coherence Mesh``,
``CreateTopology`` and ``CreateGeometry``, forms one surface loop and
one volume, and assigns Physical Surfaces ``wall``, ``inflow`` and
``outflow`` plus Physical Volume ``fluid_domain``. This is the same
geometry-independent region-aware adapter used for generic cut surfaces;
the converter does not read the pipe-bend centreline equation, bend radius or
angle.

The mesh-project importer copies the geometry manifest, region manifest and
three patch files together, preserving this contract when the generated STL
is opened from the Geometry Library. The C++17 geometry generator currently
writes an equivalent analytic STL and a parameter/patch-count report; the
region-manifest export is presently provided by the Python implementation.

2-D planar route
~~~~~~~~~~~~~~~~

For planar STL cases, the native path uses ``StlMeshReader``, welds
vertices with ``stlWeldTolerance``, extracts a single planar boundary,
and writes Gmsh geometry. A layout may provide ``transfinite_cells`` for
a structured/transfinite layout; omitting it selects automatic triangular
meshing. For external flow, the layout supplies a rectangular outer domain
and the STL loop becomes an internal obstacle boundary. For channel/internal
flow, the STL loop itself bounds the fluid region. Physical Curve names
describe patches, not inlet/outlet state values.

Gmsh mesh-size controls
-----------------------

Geometry tessellation and volume mesh density are separate controls:

* Geometry controls determine the STL approximation: sphere/spheroid
  subdivision or cylinder circumferential/axial divisions.
* Volume mesh controls determine Gmsh target characteristic length:
  ``mesh_size``, or ``mesh_size_min`` and ``mesh_size_max``.

Keep the STL bytes and geometry reconstruction settings fixed when performing
a volume-mesh refinement study. Gmsh target size is not an exact element
count. Record the Gmsh version, input STL checksum, size bounds, generated
node/cell counts and validation report for each run.

Gmsh-to-CFD handoff
-------------------

The expected pipeline is:

``STL + manifest -> Gmsh .geo -> MSH 4.1 -> CGNS -> OneFLOW grid -> CFD Case``

The C++ MSH reader in ``codes/geometry/src/GridConversion.cpp`` parses
ASCII MSH 4.1 (and the documented legacy 2.x path), maps supported
first-order element types, preserves Physical Group names and writes CGNS
through the unstructured-grid writer. The subsequent CGNS-to-OneFLOW stage
uses the existing grid reader and explicit boundary-region mapping. Physical
names are not boundary state data; inlet profiles, wall models and other
solver parameters remain separate Case configuration.

Execution and evidence boundaries
---------------------------------

Run registered producers through ``grid/stl_to_gmsh/run_cases.py`` or
the Mesh Generation / Mesh-to-CFD Validation Workbench. The runner stages
Cases under the user workspace (normally ``.oneflow/grid-development``);
do not execute a checked-in Case in place. Mesh-to-CFD Validation adds the
Gmsh-to-CGNS, CGNS-to-OneFLOW and short CFD runtime stages. A successful
producer proves mesh generation and structural checks only. A short CFD
smoke run proves pipeline execution, not physical accuracy, convergence or
grid independence.

Related references
------------------

* :doc:`stl_geometry_generation` — analytic asset definitions and generator.
* :doc:`stl_to_gmsh_2d` — planar STL route.
* :doc:`stl_to_gmsh_3d` — cube and closed-shell 3-D Cases.
* :doc:`stl_to_gmsh_sphere` — sphere-specific acceptance contract.
* :doc:`stl_to_gmsh_prolate_spheroid` — 6:1 spheroid Case.
* :doc:`stl_to_gmsh_commands` — direct volume-mesh command.


Generic STL surface cutting
---------------------------

The canonical internal-flow cylinder is also used as the first blind-input
fixture for generic surface cutting. This route intentionally does not read
its geometry manifest or use the cylinder equation. The cutter accepts only
STL triangles, two plane points, two plane normals and a geometric tolerance.

The implementation in grid/stl_import/surface_cut.py clips each triangle
against the inlet and outlet half-spaces, extracts one boundary loop on each
plane, triangulates both loops and returns a closed triangle shell. The inlet
normal points into the retained domain; the outlet normal points toward the
outlet. The generated caps are oriented outward.

The first implementation supports one simple loop per cut plane. Automatic
endpoint direction estimation, branching, multiple cross-section loops and
holes are later stages. The geometry-only regression test reads the cylinder
STL without opening its analytic manifest. Actual Gmsh volume generation,
CGNS handoff and the short OneFLOW CFD smoke run are performed manually from
the Workbench, not in CI.
