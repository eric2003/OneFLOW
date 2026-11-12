STL Unit Sphere to Gmsh 3D Tetrahedral Mesh
===============================================

Purpose
-------

This entry documents the closed-surface volume-meshing path for the unit
sphere Case. It is the first curved, closed 3-D STL example in the STL-to-Gmsh
progression. The objective is to make the geometry-to-volume-mesh operations,
their assumptions, commands and acceptance checks discoverable and
reproducible.

This is a mesh-generation contract. The later Mesh-to-CFD smoke route is
described separately; its wall-only sphere setup is not an external-flow
physical model or an accuracy benchmark.

Geometry and input STL
----------------------

The reference sphere is centered at the origin and has radius one:

.. math::

   x^2 + y^2 + z^2 = 1

The analytic reference measures are

.. math::

   V = \frac{4\pi}{3}, \qquad A = 4\pi

The checked-in ``input/model.stl`` is a triangulated, watertight icosphere
approximation with 1,280 outward-oriented triangular facets. It represents
one connected, genus-zero closed shell. Each undirected edge is used by
exactly two facets with consistent orientation, and the Euler characteristic
is two. The vertices lie on the unit sphere to the precision stored in the
STL. The STL is a piecewise-planar approximation, not an exact analytic CAD
surface.

The native C++ and Python reference Cases use byte-identical STL files. The
input shell must be closed, connected, consistently oriented, non-degenerate,
and manifold. Open surfaces, multiple disconnected shells, cavities,
self-intersections and non-manifold edges are outside this Case contract.

Geometry reconstruction and volume creation
-------------------------------------------

The shared geometric intent is:

1. Read the STL triangle facets and validate the shell topology and outward
   orientation before invoking Gmsh.
2. Merge the STL surface into the Gmsh model.
3. Classify the discrete facets with the configured surface and curve angle
   thresholds. The current Case uses 40 degrees for surface classification
   and 180 degrees for curve classification.
4. Reconstruct parametrizable surface patches from the classified facets.
5. Create one surface loop from all reconstructed patches.
6. Create one enclosed 3-D volume from that surface loop.
7. Declare the volume as the ``fluid_domain`` Physical Volume and the
   complete exterior shell as the ``wall`` Physical Surface.
8. Generate an unstructured first-order tetrahedral volume mesh and write
   ASCII MSH 4.1.

The physical names describe mesh regions and patches. They do not, by
themselves, provide inlet/outlet states or a complete solver boundary
condition. This sphere has one exterior wall group and no inlet, outlet or
farfield group.

Mesh-size refinement
--------------------

The checked-in Case layout starts with ``mesh_size = 0.25``. The runner's
refinement sweep keeps the STL and geometry-reconstruction settings fixed and
changes only the target characteristic length:

.. list-table::
   :header-rows: 1

   * - Level
     - Target characteristic length
   * - Coarse
     - 0.40
   * - Medium
     - 0.20
   * - Fine
     - 0.10

These are target sizes, not prescribed element counts. Gmsh version and
meshing algorithms can change the resulting node count, tetrahedron count
and connectivity. The acceptance contract therefore checks valid topology
and a refinement trend, not exact element totals or identical connectivity.

Implementations
---------------

The two registered implementations are:

* Native C++:
  ``grid/stl_to_gmsh/3d/sphere/automatic_tetrahedra/cpp``
  Runs through the standard OneFLOW Grid task and C++ Geometry conversion.
* Python reference:
  ``grid/stl_to_gmsh/3d/sphere/automatic_tetrahedra/python``
  Uses ``grid/stl_to_gmsh/converter.py`` and the shared STL surface
  utilities. It is an independent reference implementation, not a replacement
  for the native C++ path.

Both implementations consume the same STL and equivalent layout settings.
The expected mesh artifact is ``results/sphere_automatic_tetrahedra.msh``,
with a companion ``.stl-to-gmsh-report.json`` validation report.

Commands
--------

Build OneFLOW and make Gmsh available on ``PATH``. To generate only the
native C++ sphere mesh on Windows, stage the Case outside the repository and
run the OneFLOW Grid task:

.. code-block:: powershell

   $case = Join-Path $env:USERPROFILE ".oneflow/mesh-layout-test/grid/stl_to_gmsh/3d/sphere/automatic_tetrahedra/cpp"
   New-Item -ItemType Directory -Force (Split-Path $case) | Out-Null
   Copy-Item -Recurse -Force "grid/stl_to_gmsh/3d/sphere/automatic_tetrahedra/cpp" $case
   build/bin/Debug/OneFLOW.exe r $case

For the Python reference implementation, use a separate staging directory:

.. code-block:: powershell

   $case = Join-Path $env:USERPROFILE ".oneflow/mesh-layout-test/grid/stl_to_gmsh/3d/sphere/automatic_tetrahedra/python"
   New-Item -ItemType Directory -Force (Split-Path $case) | Out-Null
   Copy-Item -Recurse -Force "grid/stl_to_gmsh/3d/sphere/automatic_tetrahedra/python" $case
   python grid/stl_to_gmsh/run_python_case.py --case-directory $case --gmsh gmsh

The generated mesh and validation report are written under the staged Case's
``results`` directory. Never run either producer directly in the checked-in
source Case directory.

On Linux, the same procedure uses ``$HOME/.oneflow/mesh-layout-test`` as the
staging root and the built executable path, commonly
``build/bin/OneFLOW``:

.. code-block:: bash

   case="$HOME/.oneflow/mesh-layout-test/grid/stl_to_gmsh/3d/sphere/automatic_tetrahedra/cpp"
   mkdir -p "$(dirname "$case")"
   cp -a grid/stl_to_gmsh/3d/sphere/automatic_tetrahedra/cpp "$case"
   build/bin/OneFLOW r "$case"

To run the registered STL-to-Gmsh producer suite, including both sphere
implementations at coarse, medium and fine target sizes:

.. code-block:: console

   python grid/stl_to_gmsh/run_cases.py --oneflow-executable build/bin/Debug/OneFLOW.exe --gmsh gmsh

This runs the registered producer Cases, including the native and Python
sphere implementations at coarse, medium and fine target sizes. On Linux,
use the corresponding executable path, commonly
``build/bin/OneFLOW``. The runner stages every Case outside the repository
source tree; its default output root is
``%USERPROFILE%/.oneflow/mesh-layout-test/grid/stl_to_gmsh`` on Windows or
``$HOME/.oneflow/mesh-layout-test/grid/stl_to_gmsh`` on Linux.

To run one selected sphere implementation through the complete Mesh-to-CFD
smoke chain (STL -> Gmsh -> CGNS -> OneFLOW grid -> CFD runtime), use the
Workbench command-line adapter. For the native C++ Case:

.. code-block:: console

   python -m tools.oneflow_gui.mesh_validation --oneflow-executable build/bin/Debug/OneFLOW.exe --gmsh gmsh --mesh-case-id stl_to_gmsh_sphere_automatic_tetrahedra --max-steps 2

For the Python reference Case, replace the ID with
``stl_to_gmsh_python_sphere_automatic_tetrahedra``. This second command is
an end-to-end runtime smoke test, not a mesh-only command. Its default
workspace is under ``.oneflow/mesh-validation``, outside the source tree.

Mesh acceptance and evidence
----------------------------

The producer validates the input shell before meshing: non-degenerate
triangles, no duplicate facets, exactly two oppositely directed uses of every
edge, one connected component, genus-zero topology and positive signed
volume. The emitted MSH 4.1 mesh must contain first-order TETRA_4 volume
elements and TRI_3 exterior facets. The exterior-face incidence of the
tetrahedra must agree with the emitted boundary triangles. QUAD_4 and
HEXA_8 volume elements are rejected. The Physical Groups must have the
declared dimensions and names.

The Python integration tests exercise the reference producer with Gmsh when
available; the native C++ integration path is separately guarded by the
Gmsh integration-test build option. Catalog and Workbench tests verify that
both implementations are registered and that their validation graph includes
the 3-D CFD smoke stage. A linked test declaration is not itself a record of
a passing run; consult the Verification / Evidence tab and the CI revision
for execution status.

CFD smoke scope and limitations
-------------------------------

The sphere's closed exterior is entirely wall. There is no inlet or outlet.
The Mesh-to-CFD Validation Workbench deliberately reuses the laminar
flat-plate Case only for common solver parameters and overrides the staged
mesh dimension to three. This checks the mesh handoff and runtime execution
pipeline. It does not represent external flow around a sphere, and it does
not establish physical accuracy, convergence, conservation, grid
independence or a meaningful sphere drag coefficient.

This first sphere Case supports one connected, closed, consistently oriented
genus-zero STL shell and one tetrahedral volume. Multiple shells, cavities,
multi-volume assemblies, boundary-layer prisms, hybrid meshes, high-order
elements and exact CAD curvature preservation remain outside its scope.

Related entries
---------------

* :doc:`stl_to_gmsh_3d` — 3-D cube and sphere progression.
* :doc:`gmsh_interface` — MSH format and conversion semantics.
* :doc:`contract_workbench_information_architecture` — grouped Mesh
  Knowledge and shared contract navigation.
