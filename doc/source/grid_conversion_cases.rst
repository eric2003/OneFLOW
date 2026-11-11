Grid Conversion Cases
=====================

Purpose
-------

The conversion route templates are validated against small, deterministic
mesh inputs with distinct cell counts, adjacency patterns, element families,
and boundary semantics.

For the architectural ownership of the shared unstructured path (readers →
``UnstructuredMeshData`` → ``UnstructuredCgnsZoneWriter`` → CGNS Unstructured
Zone), including proves/does_not_prove and external manual references, see
the **Shared unstructured mesh conversion boundary** subsection of
:doc:`mesh_architecture`. The runnable GUI Cases describe workflows; the
small input fixtures described here provide route-level test coverage.

Conversion input fixtures are owned by their route templates. The Gmsh
two-cell topology fixture lives at
``grid/conversion/gmsh_to_cgns/grid/two_quads.msh`` and is consumed by the
automated Gmsh-to-CGNS-to-OneFLOW integration test. It is not a standalone
GUI Case. The Case catalog lists runnable route templates, while the
conversion matrix records fixture-level coverage.

Directory and metadata model
----------------------------

.. code-block:: text

   grid/
     case_catalog.json
     conversion/
       README.md
       cgns_to_oneflow/
       gmsh_to_cgns/
         grid/
           square.msh
           two_quads.msh
           cube_single_hex.msh

The responsibilities are deliberately separate:

* ``case_catalog.json`` registers runnable Cases for GUI and CLI discovery.
* A route template owns its runtime parameters and input fixtures.
* ``conversion_matrix.json`` records implemented routes and fixture-level
  conversion coverage.
* C++ reader and conversion tests stage generated outputs in temporary
  directories; checked-in input fixtures remain unchanged.
* This document explains topology, semantic contracts, and remaining gaps.

Gmsh-to-CGNS topology fixture: two_quads
---------------------------------------

The ``two_quads.msh`` fixture is a two-dimensional rectangle split into
two adjacent quadrilateral cells.

.. list-table::
   :header-rows: 1
   :widths: 42 20

   * - Quantity
     - Expected
   * - Topological dimension
     - 2
   * - Nodes
     - 6
   * - Top-dimensional cells
     - 2 QUAD_4
   * - Exterior boundary elements
     - 6 LINE_2
   * - Internal interface
     - 1 shared edge
   * - Gmsh element records
     - 9 (7 LINE_2 + 2 QUAD_4)

The cell connectivities are (1, 2, 5, 4) and (2, 3, 6, 5). Both cells share
edge (2, 5). The shared curve is present in the Gmsh input but has no
external Physical Group assignment. It must not be mistaken for an exterior
boundary condition during conversion.

The exterior boundary is split into six line elements and grouped into four
named Physical Groups: inflow, outflow, symmetry, and farfield. The bottom
and top groups each contain two curve entities. Both surface entities belong
to fluid_domain.

3D conversion reference fixture: cube_single_hex
------------------------------------------------

The automated Gmsh-to-CGNS conversion suite also includes
``grid/conversion/gmsh_to_cgns/grid/cube_single_hex.msh``, a deterministic
MSH 4.1 fixture with eight nodes, one HEXA_8 volume cell, and six QUAD_4
boundary elements. Its five Physical Groups are inflow (XMin), outflow
(XMax), symmetry (YMin), farfield (YMax), and wall (ZMin and ZMax). The
fixture is used to verify the 3D conversion boundary; it is not yet an STL
geometry reconstruction Case.

For unstructured CGNS zones, boundary elements have codimension one relative
to the cell dimension. The writer therefore uses ``EdgeCenter`` for 2D
boundary LINE_2 elements and ``FaceCenter`` for 3D boundary TRI_3 or QUAD_4
elements. BC_t PointList values reference the corresponding boundary-element
section, not volume-cell IDs. The 2D and 3D conversion tests independently
read the written CGNS files and verify the GridLocation and element references.
The 3D test additionally sends the CGNS file through CgnsFactory and checks
that OneFLOW reconstructs the expected 8-node, 1-cell, 6-boundary-face
unstructured grid.

Conversion route and execution
------------------------------

The two-quad fixture is consumed by
``GmshTwoQuadsOneFlowTest``. The test copies the MSH input to a temporary
directory and exercises the composed Gmsh -> CGNS -> OneFLOW route. The CGNS
interchange file and native OneFLOW grid are written only under that temporary
workspace; no generated files are written into the source tree.

This fixture is not a catalog entry and is not validated through
``tools/mesh_cases.py``. ``GmshMeshReaderTest`` checks its six nodes, two
quadrilateral cells, seven line elements, Physical Group assignments, and
ungrouped internal curve. ``GmshTwoQuadsOneFlowTest`` exercises the composed
conversion route and currently checks that non-empty CGNS and native-grid
artifacts are produced in a temporary workspace. It does not yet independently
read back and assert the two-cell CGNS topology. Separate CGNS round-trip tests
cover the triangle and single-hexahedron fixtures.

The current Gmsh conversion boundary is independent of SU2 parsing. The MSH
reader produces a GmshMesh value, the adapter maps supported first-order
cells and boundary patches into UnstructuredMeshData, and
UnstructuredCgnsZoneWriter constructs the CGNS Unstructured Zone. SU2
uses the same writer through a separate source-format adapter. The current
Plot3D path preserves structured blocks through GridMediator / CGNS
Structured Zones. An optional Plot3D-to-Unstructured route could expand
structured block connectivity into QUAD_4 / HEXA_8 cells and use the shared
writer, but it must explicitly handle IBLANK, block interfaces and boundary
semantics; that route is not yet implemented.


STL-to-Gmsh generation workflow
-------------------------------

The geometry-generation stage has native OneFLOW Cases under
`grid/stl_to_gmsh/square/single_quad/cpp` and `grid/stl_to_gmsh/two_square_cells/cpp`.
They run through the standard C++ Grid task, recover the planar STL boundary,
write Gmsh geometry, invoke the configured Gmsh CLI, and validate the MSH 4.1
output. Independent Python reference Cases live beside the native Cases under each
problem directory in `python/`. The existing external-flow Python preprocessing Case remains available at
`grid/stl_to_gmsh/script/run.py` with `config/external_flow.json`.

External process execution
--------------------------

Native grid conversion invokes external tools through the shared
``ONEFLOW::RunProcess`` service in ``codes/io``. It passes the executable
and each argument separately; it does not construct a shell command or rely on
``std::system``. Windows uses the wide-character process API and applies
Windows command-line argument escaping. POSIX uses ``posix_spawnp``, preserving
argument boundaries and searching the process PATH.

The process result records whether startup succeeded, the exit code, captured
standard output/error (merged into one diagnostic stream), and a startup or
wait error when applicable. The STL-to-Gmsh adapter includes this information,
the executable and working directory in its failure message. This is important
because a nonzero Gmsh exit code alone does not explain whether the failure
came from process startup, geometry parsing, meshing, or output writing.

The C++ ``ProcessRunnerTest`` covers argument boundaries, captured output,
nonzero exit codes, and missing executables. The native STL-to-Gmsh integration
tests use the same service for Gmsh discovery and invocation. When conversion
fails, the integration fixture preserves its temporary Case directory and
prints its location so that the generated ``.geo`` file and related artifacts
can be inspected rather than disappearing during test cleanup.

Optional external-tool test policy
-----------------------------------

Gmsh is an optional executable dependency. The ordinary OneFLOW test suite
must remain useful on machines where Gmsh is not installed. Configure the
test policy with ``ONEFLOW_EXTERNAL_TOOL_TESTS``:

.. code-block:: console

   cmake -S . -B build -DONEFLOW_EXTERNAL_TOOL_TESTS=AUTO
   cmake -S . -B build -DONEFLOW_EXTERNAL_TOOL_TESTS=ON
   cmake -S . -B build -DONEFLOW_EXTERNAL_TOOL_TESTS=OFF

The modes are:

* ``AUTO`` (default): register the Gmsh integration tests only when the
  ``gmsh`` executable is discoverable on PATH. No skipped tests are reported
  when Gmsh is absent.
* ``ON``: require Gmsh at configure time and fail configuration if it is
  unavailable. Use this in a dedicated external-integration environment.
* ``OFF``: do not register Gmsh integration tests, even if Gmsh is installed.

STL parsing, planar-boundary extraction, geometry-script construction,
MSH parsing, and ProcessRunner behavior remain core tests and do not require
Gmsh. Only tests that launch the actual Gmsh executable are conditional.
This setting currently controls Gmsh; additional external tools should be
added to the same policy only when their integration tests are introduced.

The regular Windows and Linux CI jobs explicitly use ``OFF`` so their core
test results do not depend on the runner image's optional software inventory.
A separate Windows Gmsh integration workflow installs a pinned Gmsh release
and configures the test suite with ``ON``, ensuring that real CLI invocation
is still exercised rather than silently skipped.

Coverage plan
-------------

The first route-owned two-cell fixture establishes a multi-cell 2-D baseline.
Extend the fixture set only when each new input adds a distinct parser,
topology, element-family, or boundary-semantic contract:

.. list-table::
   :header-rows: 1
   :widths: 24 32 44

   * - Case family
     - Candidate Cases
     - Main additional evidence
   * - Basic 2-D
     - single_triangle, single_quad
     - Element arity and minimal geometry
   * - Multi-cell 2-D
     - two_triangles, two_quads
     - Cell adjacency, shared-face identity and orientation
   * - Mixed 2-D
     - mixed_2d
     - Mixed cell sections and topology normalization
   * - Boundary semantics
     - boundary_groups
     - Multiple Physical Groups, group names and boundary mapping
   * - Basic 3-D
     - single_tetra, single_hexa
     - Volume cell and codimension-one boundary faces
   * - Multi-cell 3-D
     - two_tetrahedra
     - Shared triangular face and orientation
   * - Mixed 3-D
     - mixed_3d
     - Mixed volume cell sections and boundary semantics

These are planned examples, not claims of current support. Each should be
added only with a valid input file, an independently runnable Case, expected
topology, and a documentation page. Keep Gmsh MSH 4.1 ASCII as the initial
baseline; record high-order, periodic, partition and ghost limitations
explicitly until implemented.

The Case catalog and conversion matrix have different jobs. Do not duplicate
runtime parameter values in either registry. The Case's existing OneFLOW
parameter files remain the source of runtime configuration, while the matrix
describes route and coverage metadata.


Mesh-to-CFD validation chain
----------------------------

The STL-to-Gmsh producer Cases can be composed with the shared Mesh-to-CFD
Validation Workbench. The Workbench accepts a registered native C++ or
Python producer and constructs an artifact-bound graph:

.. code-block:: text

   selected STL-to-Gmsh Case
       |
       v
   Gmsh MSH 4.1 (validated)
       |
       v
   grid/conversion/gmsh_to_cgns
       |
       v
   CGNS interchange artifact
       |
       v
   grid/conversion/cgns_to_oneflow
       |
       v
   native OneFLOW grid + boundary-region map + grid-info JSON
       |
       v
   test/cfd/plateuns2dslau2 (runtime-overridden CFD smoke test)
       |
       v
   residual + restart flow + visualization

The reusable `grid/conversion/cgns_to_oneflow` Case is a conversion-stage template. The
graph runner injects the upstream CGNS artifact and binds the selected mesh
stem to the native-grid, boundary-region-map and grid-info outputs. The
existing `test/cfd/plateuns2dslau2` Case is staged as the CFD consumer; its
smoke-test iteration count, output cadence and initialization mode are
overridden only in the staged copy. This avoids maintaining a second solver
parameter tree solely for the STL workflow.

Run a selected producer from the repository root after building OneFLOW and
making Gmsh available on PATH:

.. code-block:: console

   python -m tools.oneflow_gui.mesh_validation --oneflow-executable build/bin/Debug/OneFLOW.exe --gmsh gmsh --mesh-case-id stl_to_gmsh_square_2x2 --max-steps 2

The default workspace is under
`%USERPROFILE%/.oneflow/mesh-validation/` on Windows and
`$HOME/.oneflow/mesh-validation/` on Linux. Each selected mesh Case has an
independent workspace containing the staged producer, conversion stages,
CFD consumer, logs and `mesh-validation-report.json`. Re-running one Case
recreates only its own workspace; sibling Case evidence is preserved. The
source repository is not used as a runtime workspace.

The former `grid/stl_to_gmsh/pipeline/` tree, fixed
`stl_square_to_cfd.json` workflow and `run_pipeline.py` entry point have
been retired. The standalone STL development runner now executes only the
STL-to-Gmsh producer Cases; end-to-end validation has one owner.

This is a runtime smoke test, not a numerical-accuracy proof. A successful
CFD run does not establish convergence, conservation, boundary-condition
correctness, or grid independence. Those require residual-history,
finite/positive flow-state, flux-balance, and analytical or
manufactured-solution verification.


STL cube single-hexahedron validation
-------------------------------------

The `grid/stl_to_gmsh/cube_single_hex` Case pair adds a constrained 3-D
producer to the STL-to-Gmsh collection. The native C++ and independent Python
reference implementations accept an axis-aligned cube STL with twelve
triangles and reconstruct six planar faces plus one volume. Gmsh must generate
exactly eight nodes, six QUAD_4 boundary elements and one HEXA_8 volume cell.
The boundary Physical Groups are XMin/inflow, XMax/outflow, YMin/symmetry,
YMax/farfield, and ZMin+ZMax/wall; fluid_domain is the volume group.

The Case appears in both Run Development Cases and Python Reference Cases,
and is discoverable by Mesh-to-CFD Validation through the Case catalog. Its
validation graph ends after STL -> Gmsh -> CGNS -> OneFLOW, where the
FaceCenter boundary locations and boundary-region map are checked. It does
not run the unrelated flat-plate CFD template. This prototype is limited to
axis-aligned cubes with equal edge lengths; it is not a general STL-to-hex
mesher. The unit-test suite is stored under tests/python/grid rather than in
the grid/stl_to_gmsh Case tree, and all executable Cases are staged outside
the source repository before running.


The `square/automatic_triangles` C++ and Python reference Cases are separate
catalog entries. They are discoverable in the development-case selectors and
Mesh-to-CFD Validation. Unlike `two_triangles`, the automatic case does not
promise an exact cell count: `mesh_size` controls the target refinement and
Gmsh determines the resulting unstructured TRI_3 topology.


Case construction as a reusable engineering method
---------------------------------------------------

A Development Case is more than a test mesh. It is a small, reproducible
engineering procedure with an explicit question, controlled inputs, a
construction method, a named artifact contract, and evidence at each
workflow boundary. The companion :doc:`development_case_methodology`
describes the Case template and the responsibilities of native C++ and
Python reference implementations. The GUI organizes these Cases by spatial
dimension and family; dimension is catalog metadata, not inferred from a
display title.


STL-to-Gmsh 3D cube Cases
-------------------------

The STL volume-meshing Cases are distinct from the checked-in Gmsh-to-CGNS
conversion fixture above. They share one closed cube STL and progressively
change the volume topology: one HEXA_8 cell, a fixed six-TETRA_4 subdivision,
and an automatic Gmsh tetrahedral mesh. See :doc:`stl_to_gmsh_3d` for the
input geometry, construction strategy, and acceptance contract.
