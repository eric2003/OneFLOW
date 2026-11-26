Mesh-to-CFD Validation
======================

Mesh-to-CFD Validation is available from the Workbench home screen, the Mesh
Generation hub, and the OneFLOW Tests dialog. It is separate from **Run
Development Cases**, which continues to execute individual preconfigured
Cases. The STL-to-Gmsh project page also opens this workbench for its current
MSH file; this project-scoped route does not regenerate the mesh.

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


Generic imported 3-D STL projects
---------------------------------

A mesh project imported from an arbitrary STL uses the geometry-neutral
``stl_to_gmsh_generic_3d_mesh_handoff`` validation entry. The entry is a 3-D
mesh input profile, not a handoff-only route: after MSH4 inspection,
Gmsh-to-CGNS conversion and CGNS-to-OneFLOW native-grid ingestion, it continues
into ``oneflow_cfd_solve``. The existing ``laminar_flat_plate`` CFD control
Case is reused as the generic solver smoke setup, with ``dimension`` overridden
to ``3`` in the isolated runtime Case. When the project carries a
checksum-verified ``region_manifest.json``, its named regions become required
MSH Physical Groups for that run.

The generated 90-degree pipe-bend Case shares the same geometry-neutral
region-aware converter, but its own lifecycle contract currently selects a
shorter milestone ending at CGNS conversion. The declared lifecycle scope
determines which workflow stages are run and which stages are expected to
remain ``NOT_RUN``; the generic four-stage route must not be inferred for
every geometry Case.

90-degree pipe bend
-------------------

The parameterized 90-degree pipe bend is registered as a generated
Python-reference Case in the 3-D **Mesh-to-CFD Validation** tree. The producer
runs the parameterized geometry generator and the shared region-aware Gmsh
volume mesher, then copies the resulting MSH into an isolated validation run.
Its current lifecycle scope includes geometry, named surface regions, the
tetrahedral volume mesh and Gmsh-to-CGNS conversion. The four required MSH
Physical Groups are dimension-2 ``wall``, ``inflow`` and ``outflow``,
plus dimension-3 ``fluid_domain``.

The declared milestone ends after CGNS conversion. CGNS-to-OneFLOW ingestion,
the native solver smoke test and physical CFD validation remain separate
future milestones; the current milestone is not considered failed because
those later stages are ``planned``. Mesh-quality metrics are diagnostic and
do not gate conversion. This milestone verifies software geometry/mesh/format
handoff only, not pressure loss, secondary-flow accuracy or mesh independence.

When a generated mesh project is opened directly, the project-scoped
**Mesh-to-CFD Validation...** action validates that exact project. Its workflow
scope remains governed by the selected Case lifecycle contract.

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

Its MSH check requires the two-dimensional Physical Surfaces ``wall``,
``inflow`` and ``outflow``, plus the three-dimensional Physical Volume
``fluid_domain``. Missing groups or groups declared with the wrong dimension
fail the mesh-generation stage before conversion. The validation then runs
the same four-stage graph used by other 3-D mesh Cases, including
``oneflow_cfd_solve``. The existing laminar flat-plate solver controls are
reused, but the runtime ``dimension.txt`` is overridden to ``3``; the native
grid and boundary-region map are passed into the solver stage in the
project-local validation workspace. This is deliberately a short solver
smoke test, not a claim that the flat-plate physical model represents a
validated internal-flow experiment. Mesh quality and physical accuracy remain
``NOT_CHECKED``. The Windows CI workflow also exercises the generic curved
two-loop STL through the same CFD smoke stage and keeps all run data under
the runner temporary directory.

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

Geometry lineage
----------------

Each run report includes a ``geometry_lineage`` record. When a companion
``.geometry.json`` is present, its asset ID, definition hash and generator
metadata are accepted only if the manifest schema is recognized and its STL
SHA-256 matches the actual source bytes. The report records the source STL
checksum and the derived MSH checksum, with the MSH linked to the source
representation. When produced, the CGNS conversion and OneFLOW native-grid
artifacts are also recorded with their own checksums and parent links. A
missing or mismatched manifest is not silently trusted.

Historical cube STL Cases currently have no verified original source or
generator record. Their Case Catalog entries explicitly declare provenance
``unknown``; the validation result still records the actual checked-in STL
path and checksum and the produced MSH path and checksum. A checksum proves
which bytes were used, not who authored the geometry or which tool originally
created it. The GUI shows these lineage rows alongside the workflow stages.

Persistent geometry lineage sidecars
-------------------------------------

Each generated CGNS volume mesh is accompanied by
``<mesh.cgns>.geometry-lineage.json``. It records the geometry provenance,
verified source STL (when available), derived Gmsh MSH, and a SHA-256 checksum
of the CGNS file itself. The following CGNS-to-OneFLOW stage verifies that
checksum before consuming the record and writes the inherited lineage beside
the native grid as ``<mesh.oneflow.grid>.geometry-lineage.json``. The native
grid sidecar links back to the CGNS filename and checksum.

These JSON files are companion metadata; they do not modify the CGNS or
native-grid binary formats. A missing, malformed, or checksum-mismatched
upstream record causes the lineage handoff stage to fail rather than silently
claiming an unverified chain. The geometry source may still be explicitly
``unknown``; persistence does not manufacture historical provenance.

Artifact inventory and report recovery
---------------------------------------

The run report also contains ``artifacts.files``, a deterministic inventory of
files that actually existed when the report was finalized. Each entry stores a
run-relative POSIX path, stage, display name, category, byte size and SHA-256 digest.
The report itself is excluded from this list to avoid a self-referential
checksum. Temporary atomic-write files and symlinks are not followed.

The inventory categories are Geometry, Volume Mesh, CFD Results, Reports,
Logs and Other. The GUI groups the entries by category and checks whether each
recorded file is present when the report is opened. A missing file remains
visible as Missing; the UI does not silently discard it. Paths are relative to
the report directory, so moving the complete run directory preserves the
artifact links. The console-log reference is also run-relative when the log
is stored inside the run directory. When a report is reopened, the GUI
re-bases its run directory on the report's current location, so moving the
complete run directory does not leave stale absolute paths in the result view.
The worker refreshes the inventory after writing the durable console log,
ensuring that the final report includes that log as well.

Reports created before ``artifacts.files`` was introduced remain readable:
their existing artifact fields and geometry lineage continue to be displayed.
The raw JSON report remains available for complete machine-readable evidence.


Recovering a run when its JSON report is missing
-------------------------------------------------

Each run now writes a small ``mesh-validation-run.json`` metadata record
before executing its stages. It preserves the Case ID, implementation, input
mode, CFD template, step limit, source paths, run ID, start time and planned
stage directories. This makes interrupted runs identifiable even when no
final report was written.

Use **Recover Folder...** and select the run directory when
``mesh-validation-report.json`` is missing or unreadable. Recovery reads
the per-stage ``results/mesh-validation.log`` records, the optional
``stl_to_gmsh/results/mesh-validation.json`` evidence, the durable
``mesh-validation-console.log``, and checksum-bound geometry lineage
sidecars. For older runs whose stage logs predate terminal-status recording,
the durable console log is used as a fallback source for the exact
``[stage] succeeded/failed/cancelled`` records.

A stage is marked PASS, FAIL or NOT RUN only when its own log contains the
terminal status line written by the stage runner. A mesh file, CGNS file, CFD
result or partial log alone never proves that its producing stage completed
successfully; such stages remain NOT CHECKED. Verified geometry sidecars are
accepted only when their artifact checksum still matches.

The recovered view is explicitly labelled RECOVERED. It does not invent the
original overall exit code, run completion time, missing stage outcomes, or
physical-accuracy verdict. The reconstructed report is an in-memory evidence
view; recovery does not overwrite or recreate the missing original report.


Run history and retention
-------------------------

Each new validation execution is stored in its own timestamped run directory below the Case validation workspace. Starting a new run no longer deletes previous reports, logs, meshes or solver outputs. The **Run History** tab lists completed reports and interrupted runs with readable run metadata; double-clicking an entry opens the recorded report or reconstructs an evidence-only view from its durable stage logs.

The history list is a filesystem view, not a second source of truth. It is rebuilt by scanning run directories, so deleting or corrupting an optional index cannot orphan a run. A completed report determines its recorded status; a directory with metadata or stage logs but no readable report is labelled ``RECOVERABLE`` until opened and reconstructed. Older runs stored directly in the validation workspace remain discoverable. Run directories are retained until the user explicitly removes them outside the application; automatic cleanup is intentionally not performed.
