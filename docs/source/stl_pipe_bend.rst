90-Degree Pipe-Bend STL Geometry
=================================

Purpose and status
------------------

The first curved-pipe geometry asset is a parameterized 90-degree circular
pipe bend under
``grid/stl_to_gmsh/3d/pipe_bend/90_degree/geometry``. It extends the
existing sphere, spheroid and straight-cylinder geometry library.

The Python geometry generator now exports a closed STL plus three separate
boundary-patch STL files and a region manifest. These feed the existing
generic region-aware STL-to-Gmsh route; the converter does not know that the
input is a pipe bend. The automatic-tetrahedra producer is registered as a runnable Case in the
3-D Mesh-to-CFD Validation tree under **STL to Gmsh**. It generates the STL,
region manifest and Gmsh MSH4 before continuing through CGNS, OneFLOW and the
configured CFD smoke test. Physical-accuracy evidence remains outstanding.

Canonical parameter set
-----------------------

The default pipe radius is 0.5 m (diameter 1 m), the centreline bend radius
is 1.0 m, the bend angle is 90 degrees, and the inlet/outlet straight
extensions are each 1.5 m. Thus

.. math::

   R_c/D = 2, \qquad L_{in}/D=L_{out}/D=3.

The surface discretization defaults are 48 circumferential divisions,
24 bend intervals, and 12 intervals on each straight section. These are
surface approximation controls, independent of the later Gmsh volume mesh
characteristic length.

Centreline and swept surface
----------------------------

Take the inlet direction as +X and the bend plane as X-Y. For the circular
arc, parameterized by :math:`0 \leq \phi \leq \theta`,

.. math::

   \boldsymbol C(\phi)=
   \begin{pmatrix}
     R_c\sin\phi \\
     R_c(1-\cos\phi) \\
     0
   \end{pmatrix}.

The unit tangent, in-plane normal and fixed binormal are

.. math::

   \boldsymbol T=(\cos\phi,\sin\phi,0),\quad
   \boldsymbol N=(-\sin\phi,\cos\phi,0),\quad
   \boldsymbol B=(0,0,1).

The circular cross-section is transported along this frame:

.. math::

   \boldsymbol P(\phi,\alpha)=
   \boldsymbol C(\phi)+r\cos\alpha\,\boldsymbol N(\phi)
   +r\sin\alpha\,\boldsymbol B.

The inlet straight section is tangent to +X at the start of the arc. The
outlet straight section begins at the arc endpoint and follows its endpoint
tangent. The frame is continuous at both joins. This is a swept circular
section along a piecewise straight/circular centreline, not a Boolean union
of a straight cylinder and a torus.

Analytic volume reference
------------------------

For a constant circular section and a non-self-intersecting bend with
:math:`R_c>r`, the enclosed fluid volume is

.. math::

   V_{analytic}=\pi r^2
   \left(L_{in}+L_{out}+R_c\theta\right),

where :math:`\theta` is in radians. The STL surface is piecewise planar,
so its signed volume is compared with this analytic reference using a
relative tolerance. The reference is a geometry check, not a CFD result.

Surface patches and boundary meaning
------------------------------------

The closed shell contains three distinct patches:

* ``wall``: the swept lateral surface;
* ``inflow``: the inlet disk;
* ``outflow``: the outlet disk.

The two disks close the surface for volume reconstruction. They are not
solid walls. The intended Gmsh contract is dimension-2 Physical Surfaces
``wall``, ``inflow``, ``outflow`` and one dimension-3
Physical Volume ``fluid_domain``. STL does not store Gmsh Physical
Groups, so patch identity must travel in a checksum-bound sidecar manifest.
The current importer discovers ``region_manifest.json`` beside the STL; keep
one region-aware generated asset per output directory. The GUI's default
output directory is already dedicated to this geometry.

Python and C++ implementations
------------------------------

Both generators write a closed STL, a checksum-bound
``oneflow.stl-geometry-manifest`` sidecar, three named patch STL files
(``wall``, ``inflow``, ``outflow``) and a checksum-bound
``region_manifest.json`` using the shared
``oneflow.stl-surface-regions`` contract. The Python generator writes
binary STL; the C++17 generator writes ASCII STL. Their STL bytes are not
expected to match, but both record the same geometric parameters, units,
surface counts and intended Gmsh Physical Groups. The C++ manifest also
records the generator source hash.

The generic converter verifies that the patch triangles form an exact
geometric partition of the closed STL, not merely a matching triangle-count
sum. It consumes only the STL patches and region manifest: it merges the
patches, coheres shared nodes, creates discrete geometry and forms the closed
fluid volume without consulting the bend equations. Both geometry generators
therefore use the same downstream meshing path. The GUI geometry editor
currently invokes the Python generator; the generated Case also exposes a
native C++17 producer that writes the same patch/manifest contract before
calling the shared Gmsh volume mesher.

Automated C++/Python parity tests verify both manifests and patch checksums,
then compare closed topology, component and vertex/triangle counts, surface
area, signed volume and coordinate bounds. A separate real-Gmsh integration
test runs the C++ producer through the shared volume mesher and checks the
four Physical Groups and non-empty TETRA_4 volume mesh in the Linux/Windows
CI matrix. These checks establish geometry and mesh-production parity, not
CGNS conversion, OneFLOW ingestion or CFD physical accuracy.

The geometry is available in **Mesh Generation -> Generate STL Geometry**
and in the **Geometry Library** parameter editor as ``90° Pipe Bend — Internal
Flow``. It is also a first-class **STL to Gmsh** Geometry preset: selecting
``90-Degree Pipe Bend STL`` generates the default checksum-bound STL and
region manifest in the user mesh-generation workspace, then assigns the
registered generated Case to the validation handoff. The volume-mesh density
remains an independent control from the surface discretization parameters. The editor exposes pipe radius, centreline bend radius, bend angle,
inlet/outlet lengths and each surface discretization control. It writes the STL
and manifest/report to the user's mesh-generation workspace. The **Open
Generated STL in STL to Gmsh** action is available: mesh-project import copies
the geometry manifest, region manifest and three patch files, so the existing
generic region-aware converter retains the inlet/outlet/wall Physical Groups.
The **Tests -> Run Development Cases** selector also offers the pipe bend as
a Python reference Case when **Python reference cases** is enabled (it remains
separate from the default native C++ selection). The standalone
**Mesh-to-CFD Validation** tree lists the pipe bend under **3D → STL to Gmsh**
as a generated Python-reference Case. Selecting it runs
the parameterized STL generator, region-aware Gmsh volume mesher, MSH4 group
checks, CGNS conversion, OneFLOW ingestion and CFD smoke test. The project
button remains available for validation from a generated mesh project. The
same validation workbench is also available from the Workbench home screen and
the Mesh Generation hub.

Both generators use the same default dimensions, discretization and declared
unit metadata. Outputs
belong under ``%USERPROFILE%/.oneflow/mesh-generation/generated_stl/``
(or ``$HOME/.oneflow/...``), never in the repository source tree.

Python command:

.. code-block:: console

   python grid/stl_to_gmsh/3d/pipe_bend/90_degree/geometry/python/generate_pipe_bend.py --output "$HOME/.oneflow/mesh-generation/generated_stl/pipe_bend_90_degree/pipe_bend_90_degree.stl"

C++ command on Windows (Visual Studio Developer PowerShell):

.. code-block:: console

   cl /std:c++17 /EHsc /O2 /Z7 /Fo"$env:USERPROFILE\.oneflow\mesh-generation\pipe_bend_generator.obj" /Fe"$env:USERPROFILE\.oneflow\mesh-generation\pipe_bend_generator.exe" grid\stl_to_gmsh\3d\pipe_bend\90_degree\geometry\cpp\generate_pipe_bend.cpp
   & "$env:USERPROFILE\.oneflow\mesh-generation\pipe_bend_generator.exe" --output "$env:USERPROFILE\.oneflow\mesh-generation\generated_stl\pipe_bend_90_degree\pipe_bend_90_degree_cpp.stl"

On Linux, use a C++17 compiler:

.. code-block:: console

   c++ -std=c++17 -O2 grid/stl_to_gmsh/3d/pipe_bend/90_degree/geometry/cpp/generate_pipe_bend.cpp -o "$HOME/.oneflow/mesh-generation/pipe_bend_generator"
   "$HOME/.oneflow/mesh-generation/pipe_bend_generator" --output "$HOME/.oneflow/mesh-generation/generated_stl/pipe_bend_90_degree/pipe_bend_90_degree_cpp.stl"

Verification and next stage
---------------------------

Run the geometry regression tests:

.. code-block:: console

   python -m unittest grid.tests.python.stl_to_gmsh.test_pipe_bend_geometry

Acceptance checks include one connected shell, watertight manifold edges,
consistent orientation, positive signed volume, relative agreement with the
analytic volume, and complete patch coverage. The default Python surface
contains 2,354 vertices and 4,704 triangles: 2,304 wall facets, 48 inlet
facets and 48 outlet facets.

A real Gmsh volume-mesh integration test is registered in the existing
Linux/Windows Gmsh 4.15.2 CI matrix. The ``automatic_tetrahedra`` mesh-producer
Case is registered in the 3-D **STL to Gmsh** list. The next verification step
is to execute the full producer-to-CFD graph and retain MSH 4.1 quality,
CGNS, native OneFLOW grid and solver smoke-test evidence. The geometry generator and the generic adapter remain
independent: changes to bend angle or centreline radius must not require any
change to the converter. A short solver smoke test must remain distinct from
physical validation of secondary flow in a curved pipe.

Related references
------------------

* :doc:`stl_geometry_generation` — existing sphere, spheroid and straight
  cylinder definitions.
* :doc:`mesh_stl_geometry_pipeline` — STL validation, Gmsh reconstruction
  and CFD handoff architecture.
* ``grid/stl_to_gmsh/3d/pipe_bend/90_degree/geometry/README.md`` —
  canonical parameters and run commands.
