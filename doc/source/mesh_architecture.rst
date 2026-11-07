OneFLOW Mesh Architecture
=========================

Purpose and status
------------------

This page is the curated overview of OneFLOW's mesh-data architecture. It brings together the repository's overlapping mesh design notes while keeping implementation facts separate from proposed interfaces.

The architecture is being developed incrementally. A proposed export level or adapter is not evidence that the corresponding implementation is complete. Current code contracts and tests take precedence over target designs.

Design goals
------------

* Inspect a mesh without starting the full solver or a graphical workbench.
* Extract a small, independently testable summary from a live runtime mesh.
* Separate mesh semantics from JSON, VTK, CGNS, Tecplot, or GUI-specific representations.
* Preserve geometry, topology, boundary meaning, and interface relationships when a format claims to support them.
* Test conversion through semantic round trips, not merely successful file creation.
* Keep derived connectivity and geometric caches distinct from authoritative source data, with explicit rebuild rules.

Information model
-----------------

The following levels provide one consistent vocabulary for the design. They describe information content, not a promise that every level is already exported.

.. list-table::
   :header-rows: 1
   :widths: 12 28 38 22

   * - Level
     - Content
     - Typical use
     - Guarantee
   * - L0 Summary
     - Dimension, counts, element-type counts, coordinate bounds, schema/version metadata
     - Quick inspection, diagnostics, inventory, UI previews
     - Inspection only; cannot reconstruct a mesh
   * - L1 Geometry and cell topology
     - Node coordinates, cell types, ordered cell-to-node connectivity where applicable, identifier conventions
     - Basic geometric reconstruction and cell-oriented algorithms
     - Reconstructs only the declared geometry and cell topology
   * - L2 Face topology and boundary semantics
     - Ordered face-to-node connectivity, owner/neighbour cells, boundary patches and physical labels, orientation conventions
     - Face-based finite-volume traversal, boundary checks, flux-compatible reconstruction
     - Preserves these semantics only when all required fields are present
   * - L3 Partition and interfaces
     - Partition/zone identity, local-to-global maps, interface-face correspondence, ghost/remote-neighbour metadata
     - Distributed mesh reconstruction and interface validation
     - Requires complete, consistently scoped maps
   * - L4 Derived data and fields
     - Centers, normals, measures, adjacency caches, search structures, solution fields and residuals
     - Visualization, diagnostics, or restoring selected solver output
     - Each item must state whether it is canonical input, reproducible derived data, or optional output

Provenance such as generator parameters, source fingerprints, conversion options, and software versions is cross-cutting metadata. It should be recorded when needed for reproducibility rather than confused with mesh topology.

Canonical data and derived data
-------------------------------

Prefer one authoritative representation for each fact.

* Node coordinates and source connectivity are canonical mesh data.
* In a face-based finite-volume model, ordered face-node connectivity and owner/neighbour relationships carry orientation and traversal meaning.
* Boundary regions and physical labels are semantic data, not display colors.
* Cell-to-face and cell-to-cell adjacency, centers, normals, areas, volumes, and search structures are derived or cached unless a specific format contract says otherwise.
* Ghost entities and partition interfaces belong to execution/distribution semantics; they must not be confused with physical cells in the source mesh.
* Generator settings describe provenance, not the generated mesh itself.

The runtime model may expose both cell-oriented and face-oriented data for algorithmic reasons. Export code should adapt from a stable read-only view or value snapshot instead of making JSON or visualization libraries part of core mesh classes.

Snapshot and adapter boundary
-----------------------------

The intended flow is:

.. code-block:: text

   Runtime Grid / UnsGrid
       -> read-only mesh view
       -> validated value snapshot
       -> format adapters
            -> versioned diagnostic package
            -> VTK XML / visualization output
            -> CGNS or other interchange formats
            -> Python inspection and visualization

L0 summary output is the first practical step. It is metadata only: it does not contain the coordinates and connectivity required to reconstruct a mesh. A richer package should declare its schema version, available levels, ID base, ID scope, units when known, and unsupported or lossy fields.

JSON is suitable for manifests and small metadata. Large numeric arrays may need compact typed payloads. JSON-library types should remain inside serialization adapters rather than leaking into the canonical mesh model.

MeshSnapshot v1 interchange contract
--------------------------------------

The first executable interchange contract is available in
`tools/mesh_contracts/snapshot.py`. It is deliberately a tooling-side
manifest validator, not a replacement for the production C++ mesh container
or a claim that runtime export is complete.

The version-1 neutral data shape contains:

* dimensioned nodes with positive IDs and finite coordinates;
* cells with a declared element type and ordered node IDs;
* faces with ordered node IDs, an owner cell, and an optional neighbour cell;
* named region groups over cells and boundary/interface groups over faces;
* source-format provenance.

The validator accepts sparse source IDs: identity is not inferred from array
position. It rejects duplicate IDs, invalid coordinate dimensions, non-finite
coordinates, unsupported element arity/dimension, dangling connectivity,
invalid face ownership, and group membership that crosses entity kinds.
Connectivity order is retained; it is not sorted or normalized.

The current contract does **not** establish geometric validity, cell closure,
face orientation against the solver convention, complete cell-to-face
incidence, boundary-condition meaning, partition/ghost semantics, or
round-trip fidelity. A Gmsh Physical Group name is not itself a OneFLOW
boundary condition. The adapter must preserve its name, dimension and
membership first; a separate explicit mapping contract must assign solver
semantics. The current manifest also does not carry source CAD entity identity
or arbitrary high-order element interpolation nodes.

Run the focused contract tests from the repository root:

.. code-block:: bash

   python -m unittest tools.mesh_contracts.test_snapshot -v

The same tests are included in Linux and Windows Python CI. They validate the
manifest contract only; they are not an experimental Gmsh Case regression.
The latter remains a GUI manual test and is intentionally excluded from
automatic GitHub execution.

The intended adapter flow is:

.. code-block:: text

   SU2 reader -----------+
                         |
   Gmsh MSH reader ------+--> UnstructuredMeshData
                         |          |
   Other unstructured --+          v
      format readers          UnstructuredCgnsZoneWriter
                                    |
                                    v
                            CGNS Unstructured Zone

   Plot3D reader -> structured GridMediator / block model
                              |
                              v
                       CGNS Structured Zone

   CGNS reader / other importers -> OneFLOW Mesh Core


A future production C++ snapshot should implement the same semantics without
making the Python validator or JSON representation a runtime dependency.
The first C++ adapter should be read-only and value-based; it should not expose
Gmsh handles or let format-specific IDs leak into solver-owned identity.

Mesh Contract Workbench
------------------------

The GUI exposes Mesh as a domain independent of Formula Knowledge:
`Contract Workbench -> Mesh`. It uses the shared information architecture
(Overview, Principle, Reference, Architecture, Implementation, Code,
Documentation, Verification / Evidence) while retaining mesh-specific
topology content.

The Topology tab displays Nodes, Cells, Faces, and Groups from a loaded
MeshSnapshot JSON artifact. It is an inspection view, not a topology editor.

The Code tab now previews the current Python contract validator and the C++
Gmsh conversion implementation directly from the repository. The preview is
bounded for large source files; the full paths remain available for opening in
the configured source editor. The Documentation tab displays the contents of
this architecture page and the Gmsh interface page rather than showing only
their filenames. These are read-only source-backed views, not duplicated
Workbench-owned documentation.

The current Gmsh C++ reader still validates and summarizes MSH 4.1 records;
it does not yet export the parsed node coordinates and element connectivity
into MeshSnapshot. Therefore, loading the bundled JSON fixture demonstrates
the topology browser, but it is not evidence that the Gmsh parser has populated
that fixture. The next adapter milestone is a C++ value-based snapshot
exporter, followed by a GUI action that consumes its generated artifact.
The bundled `tools/mesh_contracts/examples/square_snapshot.json` is a small
manual fixture for checking the UI.

Evidence is intentionally split into two independent panels:

* **MeshSnapshot contract validation** reports whether the loaded manifest
  satisfies the currently implemented schema checks.
* **Format conversion execution** remains NOT RUN unless a real converter
  execution record is supplied. Opening or validating JSON never counts as
  running Gmsh, CGNS, or a OneFLOW conversion.

The current Python validator is a first-stage tooling check. It is not a
second production mesh implementation. The target remains a read-only,
value-based C++ MeshSnapshot exported from the OneFLOW runtime mesh, with
format-specific adapters on the boundary. No Python data structure should
become a solver runtime dependency.

Spatial search and point-location boundary
--------------------------------------------

Point location is a geometry service, not an ADT-specific service. `PointLocator`
combines several responsibilities, while its candidate-search implementation is now
hidden behind a private interface in `PointLocator.cpp`:

* point-coordinate storage and stable point IDs;
* duplicate detection and the tolerance policy;
* broad-phase candidate retrieval through `HXAdtTree`;
* narrow-phase Euclidean distance checks and nearest-match selection.

These responsibilities are related but should not all depend on the same concrete
implementation. The ADT is an implementation detail of candidate retrieval. A future
spatial-index contract should expose only the operations needed by geometry code, such as
rebuilding an index, inserting a point with a payload, and querying candidates in an
axis-aligned region or within a radius. Tree nodes, ownership handles, and third-party
types should not cross that boundary.

The implemented dependency direction is:

.. code-block:: text

   PointLocator
       -> spatial-index contract
            -> OneFLOW ADT adapter (existing baseline)
            -> nanoflann adapter (optional candidate backend)

The contract must specify more than matching function names. It must state:

* whether the index owns point data or only references a coordinate dataset;
* whether insertions after construction are supported, and how index rebuilds work;
* the dimensionality and coordinate lifetime requirements;
* whether a query returns all candidates in a region, nearest neighbours, or radius matches;
* that candidate ordering is not semantically significant unless explicitly documented;
* how IDs and payloads are preserved across rebuilds;
* how tolerance-boundary cases are treated by the final geometry predicate.

The current `PointLocator` performs a broad-phase axis-aligned box query and then a
narrow-phase squared-distance check. For tolerance `t`, the candidate window is
`[x-t, x+t] x [y-t, y+t] x [z-t, z+t]`; the final predicate accepts a candidate only
when its squared distance is strictly less than `t*t`. The broad phase may return false
positives, but it must not omit a point that satisfies the narrow-phase predicate. Keeping
the final predicate in `PointLocator` prevents backend-specific boundary semantics from
silently changing deduplication behavior.

nanoflann is worth evaluating as an alternative backend, not assuming as a drop-in
replacement. Its common use is nearest-neighbour and radius search over a dataset adaptor;
OneFLOW currently inserts points incrementally into the ADT. The adapter design therefore
needs to decide whether to use a supported dynamic index, rebuild a static index at defined
points, or retain the ADT for incremental construction. Compare real OneFLOW workloads:
incremental point insertion, duplicate-heavy interface matching, query throughput, build
cost, memory use, and result equivalence. Do not select a backend from synthetic query
speed alone.

The first migration step is implemented: `PointLocator` talks to a private index contract,
and an ADT adapter translates region queries into point IDs. The narrow-phase distance
predicate now reads canonical coordinates from the locator's coordinate store instead of
reading coordinates from ADT nodes. This keeps match policy independent of tree-node layout.

The remaining migration steps are:

#. Add contract tests that can be run against each backend adapter.
#. Implement a nanoflann adapter only after choosing and testing an incremental-update policy.
#. Add representative benchmarks and compare full matching results, including tolerance
   boundaries and duplicate insertion order.
#. Keep nanoflann optional until dependency licensing, packaging, platform builds, and
   performance evidence are reviewed.

The current production backend remains ADT-backed. The new boundary is an internal seam,
not yet a public plugin API; it allows the implementation to evolve without requiring
solver, mesh-container, or point-ID changes.

Topology invariants and rebuild rules
-------------------------------------

Mesh topology checks should cover at least:

* Every node, cell, face, and region reference is in range.
* Internal faces reference the expected adjacent real cells.
* Boundary faces have explicit boundary classification and a well-defined transition to ghost-cell representation.
* Face-node ordering and owner/neighbour ordering obey the solver's orientation convention.
* Face-indexed metadata remains aligned when faces are reordered.
* Interface and partition mappings use documented identifier scopes.
* Round-trip conversion preserves coordinates, cell types, ordered connectivity, face adjacency, boundary assignments, and interface maps whenever the format claims to preserve them.

Derived adjacency is not automatically coherent if callers can mutate source arrays directly. Until topology mutations are centralized or versioned comprehensively, callers must explicitly rebuild dependent data after changing source connectivity. Documentation must not imply automatic cache invalidation where the API does not enforce it.

Validation strategy
-------------------

Use small hand-constructed meshes with known answers before relying on full CFD cases.

#. Unit tests for summary extraction and schema/version handling.
#. Geometry tests independent of solver startup and global configuration.
#. Topology tests for connectivity ranges, face ownership, orientation, and boundary/interior distinctions.
#. Export tests for deterministic ordering, ID-base conversion, and actionable errors.
#. Round-trip tests that compare semantic invariants rather than byte-for-byte file equality.
#. Integration checks with an independent consumer when practical.
#. CFD regressions for end-to-end numerical behavior; these complement rather than replace local mesh tests.

Implementation boundary
-----------------------

The repository currently distinguishes a lightweight `GridSummary` contract and its JSON writer from runtime mesh-info output. Both describe metadata rather than a complete, reconstructable mesh. The next steps are to keep the summary contract independent, clarify the runtime-to-snapshot boundary, then add geometry, topology, boundary, and interface views with focused tests.

Do not replace existing mesh containers or all exporters in one pass. Evolve one tested boundary at a time.

Related engineering records
---------------------------

Detailed design and audit notes remain under the documentation engineering area and should be reviewed alongside code changes:

* `doc/architecture/`: dated architecture audits and stable design notes.
* `doc/plans/oneflow-documentation-consolidation-plan.md`: consolidation inventory and next steps.
* Root-level `GRID_TOPOLOGY.md`, `MESH_DATA_AND_VISUALIZATION.md`, and the older top-level `docs/` mesh notes remain source material for the content review; they are not separate destinations for new documentation.


Mesh Knowledge catalog
----------------------

Mesh in Contract Workbench is a grouped knowledge catalog rather than one
static conversion page. Versioned entries live under
`tools/verification/contracts/mesh/` and use the shared contract envelope:
overview, principle, references, architecture, implementation, code,
documentation, claims, limitations and verification evidence.

The initial catalog separates **Mesh Core Concepts** from **Interchange
Formats**. Core concepts cover identity, coordinates, ordered connectivity,
faces, adjacency, boundaries, regions, physical groups, units, provenance
and quality. Format entries cover Gmsh MSH, Plot3D structured grids and CGNS.
A separate topology inspector remains available for inspecting a loaded
MeshSnapshot fixture; it is not the primary Mesh knowledge interface.

The dependency direction is:

.. code-block:: text

   Gmsh / Plot3D / CGNS specifications
                  |
          format knowledge records
                  |
        format-specific adapters
                  |
       neutral OneFLOW MeshSnapshot
                  |
           OneFLOW Mesh Core
                  |
       solver / visualization / writers

The Mesh Core vocabulary is the stable center. Formats are boundary contracts
with different capabilities. A structured Plot3D block is not an unstructured
Gmsh element block; CGNS zones, sections, boundary types and connectivity must
not be flattened into anonymous cells. Conversion can preserve coordinates
while losing high-order interpolation, units, group labels or interfaces.
Those losses must be declared and verified rather than silently discarded.

The Workbench search covers titles, descriptions, source paths and references.
Each catalog entry links source files, maintained documentation and declared
evidence. A declared test is not a recorded execution result; a format entry
marked proposed does not imply a production reader or writer exists.


Gmsh reader and conversion implementation milestone
----------------------------------------------------

The MSH 4.1 ASCII reader materializes coordinates, sparse source IDs,
geometric entity ownership, physical names, element types, and ordered
connectivity in a C++ GmshMesh value model. The conversion path maps
top-dimensional first-order cells into the shared ``UnstructuredMeshData``
staging model and writes the result through ``UnstructuredCgnsZoneWriter``.
The Gmsh adapter does not construct a ``Su2Grid`` object or call a
SU2-named writer method. The bundled square mesh is now
written as a real unstructured CGNS zone and independently reopened through
the CGNS API in a focused GTest.

The conversion test verifies the Base and Zone dimensions, node and cell
counts, ordered TRI_3 connectivity, coordinate arrays, the BAR_2 boundary
section, and the BC_t PointList. The test also checks the generated
human-readable report. This validates a real file round trip rather than a
placeholder or a successful writer call alone.


Gmsh MSH 4.1 read and conversion reports
----------------------------------------

The production C++ reader maintains an explicit source-node-tag to
zero-based node-position map while retaining each original Gmsh tag. This
makes sparse source IDs safe for connectivity conversion without confusing
identity with storage position.

After a successful CGNS write, the conversion entry point writes two
case-local reports under `results/`:

* `*.gmsh-read-report.json`: machine-readable source tag ranges, coordinate
  bounds, element-type counts, and geometry checks.
* `*.gmsh-report.md`: an English inspection report covering physical groups,
  geometric entities, element and node previews, connectivity, coordinate
  bounds, and geometry checks.

The Markdown report distinguishes Gmsh input element types from the CGNS
boundary section types produced by conversion. For example, a two-dimensional
Gmsh `LINE_2` boundary is represented as a CGNS `BAR_2` section. This
boundary-section inventory is based on codimension-one elements associated
with Physical Groups, matching the conversion path; it does not invent
solver-specific wall, inlet, or outlet semantics.

Triangle area is computed from the three-dimensional cross-product magnitude;
the signed XY projection is reported separately as an orientation diagnostic.
A negative XY orientation is not automatically treated as a read failure.
Degenerate triangles are counted using a scale-relative floating-point
threshold.

The report's `READ_SUCCESS` status means that the reader's structural checks
passed and that the declared PhysicalNames and Entities sections were parsed
when present. It does not by itself claim that CGNS conversion succeeded:
the reports are emitted only after the CGNS writer returns successfully, and
the independent CGNS read-back test verifies the generated file.


Gmsh-to-CGNS first Zone
-----------------------

The first production output path maps MSH 4.1 ASCII top-dimensional,
first-order cells into ``UnstructuredMeshData`` and calls
``UnstructuredCgnsZoneWriter::Write()`` to populate OneFLOW's CGNS
Unstructured Zone model. ``CgnsFactory::WriteCgnsFile()`` serializes the
prepared base and zone through the existing CGNS writer. The Gmsh reader
itself has no dependency on the Gmsh SDK or CGNS library.

The common writer was extracted from the former
``Su2Grid::FillSU2CgnsZone()`` implementation. SU2 keeps its source-format
reader and adapts its parsed data to the shared staging model before calling
the common writer. Reusing SU2 parsing or SU2 file semantics was never a
technical requirement for Gmsh; the earlier route was an implementation
shortcut that coupled unrelated responsibilities.

The current staging model is intentionally limited: it uses VTK element
identifiers and explicit CGNS BC type values for the first-order conversion
path. It is not yet the canonical Mesh Core model. Original source entity
identity, high-order interpolation data, periodic links, partitions, ghost
entities, richer boundary semantics and complete provenance remain future
extensions.

Standard Plot3D follows a separate native structured path. The current
implementation uses ``Plot3D::Plot3DToCgns()``, ``GridMediatorGroup``,
and ``CgnsFactory::DumpCgnsGrid()`` to preserve structured blocks and write
CGNS Structured Zones.

A separate, explicit Plot3D-to-CGNS-Unstructured route is also possible
without involving SU2. Its adapter would expand structured index connectivity
into ordered QUAD_4 cells for 2-D blocks or HEXA_8 cells for 3-D blocks, then
populate ``UnstructuredMeshData`` and call
``UnstructuredCgnsZoneWriter``. This is a topology transformation, not a
native unstructured Plot3D read. It must define i-j-k node ordering, degenerate
dimensions, IBLANK handling, multi-block interface policy, and the source of
boundary patch/BC semantics. Coordinates alone must not be used to invent
inflow, wall, symmetry, or other solver boundary conditions. This optional
route is not yet implemented.

For the bundled square sample, the semantic result is four nodes, two TRI_3
cells, four BAR_2 boundary elements, and ordered volume connectivity
[1, 2, 3] and [1, 3, 4]. The four lower-dimensional LINE_2 elements are
grouped through their curve-entity physical tags and written as a CGNS
boundary section. The matching BC_t point list references the boundary
element IDs after the volume-cell range. Its name is preserved from
PhysicalNames and its type is BCTypeUserDefined; this is a label, not a
solver boundary-condition assignment.

The CGNS output path keeps connectivity one-based until CGNS serialization.
The existing SU2-to-OneFLOW runtime path still requests zero-based internal
normalization. The SU2 adapter forwards this option to
UnstructuredCgnsZoneWriter::Write(). The Gmsh file writer disables
normalization so CGNS connectivity and BC element IDs remain one-based in the
serialized file.


Gmsh Physical Group boundary semantics
--------------------------------------

The bundled two-dimensional Gmsh conversion Case now uses four named curve
Physical Groups: `inflow` on the left edge, `outflow` on the right edge,
`symmetry` on the bottom edge, and `farfield` on the top edge. The
conversion adapter maps these explicit names to CGNS `BCInflow`,
`BCOutflow`, `BCSymmetryPlane`, and `BCFarfield`, respectively. The canonical input spelling is `symmetry`; misspellings are not treated as
standard boundary-condition aliases.

This is a name-based semantic mapping, not geometric inference. Unrecognized
Physical Group names remain `BCTypeUserDefined`; the converter does not
invent a wall/inlet/outlet classification from edge coordinates. The
round-trip test independently reads all four BC_t records and verifies both
their CGNS type and the edge connectivity referenced by each PointList.

Small Gmsh fixtures should be organized by the capability they prove:
parser-only meshes, supported first-order 2-D and 3-D conversion meshes,
higher-order read-only meshes, and malformed/unsupported inputs. Use one
cell (or the minimum connected patch needed) per element family. Reader
tests should not depend on CGNS; conversion tests should write to a temporary
directory and independently reopen the result with the CGNS API; the manual
Case should remain the end-to-end GUI and report-generation path.

Format interoperability routes
-------------------------------

Format conversion is organized as a graph of supported routes, not as a
collection of pairwise converter implementations. Each format adapter owns
its file syntax and maps into or out of a stable OneFLOW-facing mesh
representation. Shared topology, entity identity, connectivity ordering and
boundary semantics must be validated by common contracts.

CGNS is the preferred standards-based interchange hub where it provides a
useful independent serialization boundary. It is not the OneFLOW runtime
mesh model and should not leak CGNS library types into the Mesh Core. A
composed route such as Gmsh -> CGNS -> OneFLOW should write and reopen a real
CGNS file; an in-memory shortcut is not equivalent evidence.

The current route inventory is maintained in
`grid/interoperability/conversion_matrix.json`. It records implementation
status and known limitations so that documentation and test planning do not
silently imply support for every possible format pair.

Current policy:

* Keep source-format fixtures in one format-owned directory. Route cases
  reference those fixtures instead of copying them.
* Keep route-specific configuration and generated evidence under the route.
* Reuse format adapters and shared semantic contracts; do not implement a
  dedicated converter for every pair.
* Add a route only when it unlocks a real workflow or an independent
  cross-validation path.
* A successful file write is insufficient: round-trip checks compare
  dimensions, coordinates, element families, ordered connectivity, region
  membership and boundary semantics.
* Record unsupported features and lossy mappings explicitly. Physical group
  names are not solver boundary states unless an explicit mapping contract
  says so.

The Gmsh -> OneFLOW convenience route currently composes the MSH 4.1 ASCII
reader, CGNS writer, and CGNS reader. It retains the intermediate CGNS file
under the Case `results/` directory as reproducible evidence. This
deliberately tests the same on-disk path that users can invoke as two
separate conversions; it does not create a second Gmsh-specific runtime
adapter.
