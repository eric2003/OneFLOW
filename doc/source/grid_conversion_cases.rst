Grid Conversion Cases
=====================

Purpose
-------

The conversion Case collection measures more than whether one element can be
parsed. A converter is only credible when it has been exercised against
multiple concrete meshes with different cell counts, adjacency patterns,
element families, and boundary semantics.

For the architectural ownership of the shared unstructured path (readers →
``UnstructuredMeshData`` → ``UnstructuredCgnsZoneWriter`` → CGNS Unstructured
Zone), including proves/does_not_prove and external manual references, see
the **Shared unstructured mesh conversion boundary** subsection of
:doc:`mesh_architecture`. Runnable Cases here exercise that boundary; they do
not redefine it.

Runnable Cases are named and grouped by the mesh they describe, not by a
conversion route. A directory such as grid/conversion/cases/two_quads remains
meaningful if the same input is later sent through Gmsh-to-CGNS,
Gmsh-to-OneFLOW, or another supported route. This avoids multiplying
directories with gmsh_to_cgns_ prefixes and keeps geometry coverage
independent from converter coverage.

Directory and metadata model
----------------------------

.. code-block:: text

   grid/
     case_catalog.json
     interoperability/
       conversion_matrix.json
     conversion/
       README.md
       cases/
         two_quads/
           README.md
           grid/
             two_quads.msh
           script/
             control.txt
             simu.txt
             dimension.txt
             tolerence.txt
             grid.txt

The responsibilities are deliberately separate:

* case_catalog.json registers each runnable Case for the GUI and CLI.
* A Case directory owns its input mesh, runtime parameters, expected topology,
  and human-readable explanation.
* conversion_matrix.json records which source/target routes are implemented
  and which runnable Cases exercise them.
* This document explains topology, semantic contracts, and remaining gaps.

A Case must be independently copyable and runnable. Do not make it depend on
a mesh file hidden in another Case directory. Shared low-level test fixtures
may still exist for unit tests, but they are not runtime dependencies of the
runnable Case collection.

Current runnable example: two_quads
-----------------------------------

The two_quads Case is a two-dimensional rectangle split into two adjacent
quadrilateral cells.

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

Conversion route and execution
------------------------------

The current Case runtime configuration selects gmsh -> cgns. It is registered
in grid/case_catalog.json, so it appears in the standard GUI Case Workbench.
Run it manually from the GUI after building OneFLOW with CGNS support. The GUI
uses an isolated run workspace; the checked-in Case and input MSH file must
remain unchanged.

The command

.. code-block:: console

   python tools/mesh_cases.py validate two_quads

checks the catalog entry and its runtime parameter includes. It is a
preflight check, not a conversion run. The Gmsh-to-CGNS runtime example is
intentionally excluded from CTest, Python unittest, and GitHub Actions.

The expected output path configured by the Case is
results/two_quads.cgns. Inspect the CGNS Base/Zone, coordinate arrays,
quadrilateral section, boundary sections, BC names and BC types. Also inspect
the reader/conversion reports. Process success or file existence alone is not
a complete topology acceptance criterion.

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
`grid/stl_to_gmsh/square` and `grid/stl_to_gmsh/two_square_cells`.
They run through the standard C++ Grid task, recover the planar STL boundary,
write Gmsh geometry, invoke the configured Gmsh CLI, and validate the MSH 4.1
output. Independent Python reference Cases live under
`grid/stl_to_gmsh/python/`. The earlier Python prototype under
`grid/stl/stl_to_gmsh/` remains available for compatibility.

Coverage plan
-------------

The first runnable Case establishes a multi-cell 2-D baseline. Extend the
collection with independent Cases, each carrying its own mesh, configuration,
expected counts, semantic description, and documentation:

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
