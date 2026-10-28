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
