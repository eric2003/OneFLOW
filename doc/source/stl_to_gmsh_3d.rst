3D STL-to-Gmsh Cube Volume Cases
=================================

Purpose
-------

These Cases isolate three volume-mesh strategies while keeping the input
geometry fixed. All variants use the same closed, axis-aligned cube STL with
12 surface triangles, eight geometric corners, six planar square faces, and
a closed two-manifold edge incidence.

Progression
-----------

.. list-table::
   :header-rows: 1
   :widths: 24 34 42

   * - Case
     - Expected volume mesh
     - Construction contract
   * - ``3d/cube/single_hex``
     - 1 HEXA_8, 6 QUAD_4 boundary faces
     - Transfinite curves, surfaces and volume; recombine surface and volume.
   * - ``3d/cube/six_tetrahedra``
     - 6 TETRA_4, 12 TRI_3 boundary facets
     - Transfinite curves, surfaces and volume; do not recombine.
   * - ``3d/cube/automatic_tetrahedra``
     - Gmsh-generated TETRA_4 count
     - No transfinite constraints; ``mesh_size`` controls characteristic size.

The six-tetrahedron Case is the fixed-count bridge between a single
hexahedron and automatic tetrahedral meshing. Its acceptance criteria are
exact: eight nodes, twelve triangular boundary facets, and six tetrahedral
volume cells. The automatic Case intentionally does not fix node or cell
counts. It requires non-empty tetrahedral volume cells and triangular
boundary facets, and rejects QUAD_4/HEXA_8 output.

Boundary semantics
------------------

All variants preserve six independent surface Physical Groups:
``xmin_inflow``, ``xmax_outflow``, ``ymin_wall``,
``ymax_farfield``, ``zmin_symmetry``, and ``zmax_symmetry``.
The volume Physical Group is ``fluid_domain``. The groups encode patch
identity and broad BC categories; they do not define inlet/outlet state values
or establish a physical cube-flow setup.

Implementations and evidence
----------------------------

Each strategy has a native C++ Case and an independent Python reference Case.
Both read byte-identical STL inputs and equivalent layout JSON. The native
implementation exercises the standard OneFLOW Grid task and
``GridConversion.cpp`; the Python reference exercises
``grid/stl_to_gmsh/converter.py`. Both write MSH 4.1 and validate element families,
counts (exact or minimum, depending on the Case), and Physical Groups.

Run all registered STL-to-Gmsh development Cases using the isolated runner:

.. code-block:: console

   python grid/stl_to_gmsh/develop_test/run_cases.py --oneflow-executable build/bin/Debug/OneFLOW.exe --gmsh gmsh

Generated meshes and reports are staged under the user's external
``.oneflow/mesh-layout-test/grid/stl_to_gmsh` directory. Do not run a checked-in
Case in place.

The Mesh-to-CFD Validation Workbench discovers all six 3D implementations
from ``grid/case_catalog.json`` and sends their MSH 4.1 outputs through the
shared CGNS and OneFLOW import stages. The short solver smoke test checks
runtime compatibility only; it is not a cube-flow accuracy or convergence
test.

Limitations
-----------

The three Cases remain restricted to one axis-aligned cube with six planar
faces. They do not claim support for arbitrary closed STL volumes, curved CAD
surfaces, cavities, multiple volumes, or general tetrahedral/hexahedral
meshing.
