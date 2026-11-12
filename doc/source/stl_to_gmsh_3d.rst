3D STL-to-Gmsh Volume Cases
===========================

Purpose
-------

The 3D STL-to-Gmsh cases progress from a topology-controlled cube to a curved,
closed unit-sphere surface. The cube cases isolate volume element topology;
the sphere case exercises Gmsh surface classification and reconstruction from
a triangulated STL shell before generating an unstructured tetrahedral volume.

Progression
-----------

.. list-table::
   :header-rows: 1
   :widths: 28 34 38

   * - Case
     - Expected volume mesh
     - Construction contract
   * - ``3d/cube/single_hex``
     - 1 HEXA_8, 6 QUAD_4 boundary faces
     - Axis-aligned cube, transfinite curves/surfaces/volume and recombination.
   * - ``3d/cube/six_tetrahedra``
     - 6 TETRA_4, 12 TRI_3 boundary facets
     - Same cube, exact fixed tetrahedral decomposition.
   * - ``3d/cube/automatic_tetrahedra``
     - Gmsh-generated TETRA_4 count
     - Same cube, automatic unstructured volume meshing.
   * - ``3d/sphere/automatic_tetrahedra``
     - Gmsh-generated TETRA_4 count
     - Unit sphere, classify STL facets, reconstruct surface patches and mesh the enclosed volume.

The cube family preserves one geometry and six independent boundary patches
while changing only the volume-mesh strategy. The sphere adds a curved closed
surface and one physical wall group. It does not add inlet/outlet semantics or
CFD physics.

Unit-sphere geometry contract
-----------------------------

The sphere is centered at the origin with radius 1:

.. math::

   x^2 + y^2 + z^2 = 1

Its exact reference measures are:

.. math::

   V = \frac{4\pi}{3}, \qquad A = 4\pi

The checked-in STL is a watertight icosphere approximation with 1,280 outward
oriented triangular facets. It is one connected, genus-zero closed surface:
each undirected edge is shared by exactly two facets, adjacent facets have
consistent winding, and the Euler characteristic is 2. Every STL vertex lies
on the unit sphere to the precision of the STL coordinates. The STL is a
piecewise-planar approximation, not an exact analytic sphere.

Both native C++ and Python reference Cases consume byte-identical copies of
this same STL. The Gmsh geometry script merges the STL, calls
``ClassifySurfaces`` with parametrizable patch reconstruction, creates the
geometry, collects all classified surfaces into one Surface Loop, and creates
one enclosed Volume. The Physical Groups are ``wall`` (dimension 2)
and ``fluid_domain`` (dimension 3).

Mesh-size refinement
--------------------

The runner holds the STL bytes and geometry-reconstruction settings fixed and
changes only ``mesh_size``:

.. list-table::
   :header-rows: 1

   * - Level
     - Characteristic length
   * - Coarse
     - 0.40
   * - Medium
     - 0.20
   * - Fine
     - 0.10

These values are target characteristic lengths, not exact element counts.
Gmsh version, surface parametrization and meshing algorithms can affect the
resulting counts and connectivity. The Python integration test verifies that
node and TETRA_4 counts increase as the requested size decreases. Exact
connectivity and exact counts are intentionally not part of the contract.

Surface validity and volume-mesh acceptance
-------------------------------------------

Before generating the geometry, both implementations validate the input shell:
non-degenerate triangles, no duplicate facets, exactly two oppositely directed
uses of each edge, one connected component, genus-zero topology, and outward
orientation with positive signed volume. The generated MSH 4.1 mesh must
contain first-order TETRA_4 volume cells and TRI_3 boundary facets. The
tetrahedral exterior-face incidence must match the emitted boundary triangles.
QUAD_4 and HEXA_8 elements are rejected. The two Physical Groups must have the
declared dimensions and names.

The reference volume and area are useful for later geometric error checks;
this first meshing Case does not yet claim a convergence study of the
piecewise-linear STL approximation or a CFD accuracy benchmark.

Implementations and evidence
----------------------------

The native implementation exercises the standard OneFLOW Grid task and
``GridConversion.cpp``. The Python reference uses
``grid/stl_to_gmsh/converter.py``. Both write MSH 4.1 and a validation report.
The orchestrator runs the sphere pair at all three mesh sizes in an external
staging directory:

.. code-block:: console

   python grid/stl_to_gmsh/run_cases.py --oneflow-executable build/bin/Debug/OneFLOW.exe --gmsh gmsh

The staging root defaults to the user's ``.oneflow/mesh-layout-test/grid/stl_to_gmsh``
directory. Checked-in Case directories are inputs only; never run them in place.

Limitations
-----------

The cube strategies remain restricted to the documented axis-aligned cube
fixture. The sphere pipeline currently requires one connected, closed,
consistently oriented genus-zero STL shell and generates one tetrahedral
volume. Multiple shells, cavities, non-manifold input, arbitrary multi-volume
assemblies, boundary-layer prisms and curved CAD surface fidelity studies are
not yet covered. The sphere case validates geometry conversion and volume-mesh
topology; it is not yet a physical CFD validation case.
