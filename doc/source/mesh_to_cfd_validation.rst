Mesh-to-CFD Validation
======================

Mesh-to-CFD Validation is a dedicated mode in the OneFLOW Tests dialog. It is
separate from **Run Development Cases**, which continues to execute individual
preconfigured Cases.

The initial conversion route is:

.. code-block:: text

   selected STL-to-Gmsh Case -> Gmsh MSH -> Gmsh-to-CGNS Case
       -> CGNS-to-OneFLOW Case -> selected CFD template

The mesh producer can be either the native C++ OneFLOW Grid task or an
independent Python reference Case. Both implementations publish the same
Case-specific MSH artifact contract. Native Cases declare the artifact path in
``script/grid.txt``; Python Cases declare the mesh filename in
``expected.json``. The language-specific producer is selected from the Case
Catalog, while the downstream conversion and validation graph is shared.

The first CFD consumer is the existing ``test/cfd/plateuns2dslau2`` Laminar
Flat Plate Case. Each source Case is copied into a unique directory under
``~/.oneflow/mesh-validation`` (or the equivalent user-home directory on
Windows). Runtime parameter changes are applied only to those copies. Every
producer and consumer process runs with its own staged Case directory as its
working directory. Users do not edit workflow JSON or manually copy
intermediate files.

The default solver smoke test advances two steps. A square mesh with the
flat-plate template tests conversion, native-grid ingestion, initialization,
and short time advancement; it is not a physical flat-plate validation.
Each stage writes a log, and the run writes ``mesh-validation-report.json``.
Additional CFD templates, including a turbulent flat plate, can be registered
as this mode evolves.

The 6:1 prolate spheroid is registered as a native C++ Development Case and
an independent Python reference Case. Both appear in the same 3-D Mesh-to-CFD
selection tree. Its wall-only Physical Group makes this a conversion and
runtime smoke test, not an external-flow accuracy benchmark.


3-D internal-flow cylinder
--------------------------

The generated **3D Cylinder Internal Flow** mesh is available as an
existing-mesh validation input. It is recorded in the catalog's separate
``mesh_validation_inputs`` collection, not as a runnable Development Case. The dialog resolves
``~/.oneflow/mesh-generation/cylinder_internal_flow/results/cylinder_internal_flow.msh``
(or the equivalent Windows user-data path), copies it into the project's
isolated ``validation/`` workspace, and runs the MSH4 inspection,
Gmsh-to-CGNS conversion, and CGNS-to-OneFLOW ingestion stages. It does not
rerun STL/Gmsh generation or modify the source project mesh.

This entry is deliberately ``mesh_handoff``, not a CFD smoke test. Its
MSH check requires the two-dimensional Physical Surfaces ``wall``,
``inflow`` and ``outflow``, plus the three-dimensional Physical Volume
``fluid_domain``. Missing groups or groups declared with the wrong
dimension fail the mesh-generation stage before conversion. The internal-flow
cylinder has wall, inlet and outlet patch semantics; the
available laminar flat-plate consumer is not a compatible physical setup.
A dedicated internal-flow cylinder CFD Case with defined boundary states and
acceptance criteria is required before a solver run can be reported as a
meaningful CFD validation. Mesh quality and physical accuracy remain
``NOT_CHECKED``.

Boundary type mapping
---------------------

The Gmsh-to-CGNS adapter maps common Physical Group name aliases to CGNS
``BCType_t`` values. Matching is case-insensitive and ignores punctuation;
the original Physical Group name is retained as the boundary patch name.

.. list-table::
   :header-rows: 1

   * - Gmsh Physical Group name
     - CGNS ``BCType_t``
     - OneFLOW boundary category
   * - ``Inlet``, ``Inflow``, ``VelocityInlet``, ``MassFlowInlet``, ``PressureInlet``
     - ``BCInflow``
     - ``BC::INFLOW``
   * - ``Outlet``, ``Outflow``, ``PressureOutlet``
     - ``BCOutflow``
     - ``BC::OUTFLOW``
   * - ``Wall``, ``LowerWall``, ``UpperWall``, ``BodyWall``
     - ``BCWall``
     - ``BC::SOLID_SURFACE``
   * - ``NoSlipWall``, ``ViscousWall``
     - ``BCWallViscous``
     - ``BC::SOLID_SURFACE``
   * - ``SlipWall``, ``InviscidWall``
     - ``BCWallInviscid``
     - ``BC::SOLID_SURFACE``
   * - ``Symmetry``, ``SymmetryPlane``
     - ``BCSymmetryPlane``
     - ``BC::SYMMETRY``
   * - ``Farfield``, ``Freestream``
     - ``BCFarfield``
     - ``BC::FARFIELD``
   * - unrecognized names
     - ``BCTypeUserDefined``
     - ``BC::GENERIC_2`` (existing fallback)

This mapping is deliberately a boundary-type mapping, not a complete
physical boundary-state definition. ``BCInflow`` does not specify velocity,
mass flow, temperature or turbulence quantities. ``BCOutflow`` does not
specify a target pressure. Mapping ``NoSlipWall`` to the broad OneFLOW
``SOLID_SURFACE`` category does not prove that a particular viscous wall
treatment is configured. The original name and CGNS type are retained for
diagnostics and later richer contracts.

Periodic connectivity is not inferred from a name such as ``Periodic``:
CGNS periodic interfaces are represented through zone connectivity metadata,
not as a ``BCType_t`` member in the enum. Neumann, Dirichlet, degenerate
and other CGNS types remain unmapped until their OneFLOW semantics and
numerical treatment are explicitly defined.

The Mesh-to-CFD dialog displays this common mapping as a quick reference.
It is not a compatibility approval: physical problem type, required boundary
state data and solver support must still be checked before a CFD template runs.


Run-result reporting
--------------------

Every Mesh-to-CFD run writes mesh-validation-report.json using the shared
oneflow.run-report version 1 envelope. The top-level status and exit_code
give the overall verdict; stages separates mesh generation, Gmsh-to-CGNS
conversion, CGNS-to-OneFLOW ingestion and the CFD smoke run. Blocked stages
are NOT_RUN, while mesh quality and physical CFD accuracy are NOT_CHECKED
until explicit acceptance criteria exist. The report also retains the
workflow's selected Case, input mode, run directory, per-stage messages
and compatibility metadata.

The command prints a short SUCCESS/FAILED summary with stage states, mesh
counts, output locations and the report path. A successful short solver
smoke run is not a physical-accuracy validation.


The dialog also presents a compact per-Case result tree: one row for each
workflow stage, its status, and its short diagnostic message. The full
terminal log remains available separately from that summary.
