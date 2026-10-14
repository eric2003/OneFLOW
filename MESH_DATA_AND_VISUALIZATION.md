# Mesh Data, Reconstruction, and Visualization Architecture

This document defines an incremental architecture for exporting, rebuilding, and visualizing OneFLOW meshes. It is a design contract, not a claim that every proposed adapter already exists.

## 1. Goals

- Separate mesh meaning from any one in-memory class or file format.
- Support lightweight inspection as well as lossless reconstruction.
- Make data loss explicit when exporting to formats with weaker semantics.
- Allow native viewers, third-party applications, and Python tools to consume the same mesh.
- Test architecture through round-trip reconstruction, not only through successful file creation.

## 2. Information levels

### L0: Summary

Suggested fields:

- schema version and mesh identifier
- topological dimension and coordinate dimension
- node, cell, face, boundary-face, and zone counts
- cell-type histogram
- coordinate bounds
- available data levels and export formats

A summary is for logging, dashboards, and quick validation. It is not sufficient to reconstruct a mesh.

### L1: Geometry and cell topology

Required for reconstructing the basic mesh:

- node coordinates, with explicit coordinate dimension and units
- cell type for every cell, or a lossless cell-type grouping
- ordered cell-to-node connectivity where the source representation provides it
- explicit index base and ID scope
- mesh and zone identity

For a polyhedral finite-volume mesh, cell-to-node connectivity alone may not preserve face orientation. L1 therefore does not replace L2 when finite-volume semantics matter.

### L2: Face topology and boundary semantics

Required to preserve the current face-based finite-volume model:

- ordered face-to-node connectivity
- face owner/left cell and neighbour/right cell
- explicit boundary-face classification
- boundary patch/region identity and physical boundary type
- face orientation convention
- interface classification and mappings when present

The export must distinguish a physical boundary from an internal face, a partition interface, and a ghost-cell relationship. A sentinel integer alone is not a portable format contract.

### L3: Partition and interface information

Required when reconstructing a distributed mesh layout:

- partition and zone identifiers
- local-to-global node, face, and cell mappings where available
- interface-face correspondence
- ghost-cell ownership and remote-neighbour mapping
- the scope of each identifier (global, zone-local, or partition-local)

A partitioned mesh without these maps may still be visualizable, but it cannot necessarily be reassembled into the original distributed layout.

### L4: Derived geometry and solution data

Examples include cell volumes, face centers, normals, adjacency caches, search trees, flow fields, and residuals.

Every field should declare whether it is canonical input, reproducible derived data, or optional solver output. Reproducible derived data should normally be rebuilt from canonical geometry and topology. Solution fields need association metadata such as node-centered or cell-centered, component count, units, and time/iteration.

## 3. Canonical model and adapters

The proposed flow is:

`OneFLOW grid -> read-only mesh view -> validated mesh snapshot -> format adapters`

The read-only view exposes semantic concepts rather than internal storage details. A snapshot owns a stable copy suitable for serialization and validation. Adapters translate the snapshot to and from external formats.

Do not make every format adapter depend directly on `UnsGrid`, `FaceTopo`, or other mutable implementation classes. Keep conversion code at the boundary. An adapter must report unsupported or lossy fields instead of silently dropping them.

Recommended initial adapters:

- **OneFLOW diagnostic package:** a versioned manifest plus explicit topology/geometry payloads; useful for reproducible round-trip tests.
- **VTK XML:** suitable for visualization and tools such as ParaView and PyVista. Preserve supported cell types and point/cell data; record limitations for face-based boundary and partition semantics.
- **CGNS:** suitable for CFD mesh interchange when the required zone, section, and boundary-condition concepts map cleanly.
- **Tecplot-compatible output:** useful for selected visualization workflows, with documented limitations for general polyhedral and interface data.
- **Python:** a small reader for the diagnostic package, with optional conversion to VTK/PyVista. Matplotlib can display simple 1D/2D cases; PyVista/VTK is more suitable for general 3D meshes.

These are choices for staged implementation, not a requirement to implement every format at once.

## 4. Can exported information reconstruct a mesh?

Use the following levels of guarantee:

- **Summary export:** inspection only; no reconstruction guarantee.
- **L1 export:** reconstructs geometry and the supported cell topology.
- **L1 + L2 export:** aims to reconstruct the same finite-volume face topology and boundary semantics.
- **L1 + L2 + L3 export:** aims to reconstruct partition and interface relationships as well.
- **L4 export:** optionally restores fields and cached data; derived caches may instead be rebuilt.

A round-trip test should compare semantic invariants, not byte-for-byte file equality:

1. node, cell, face, zone, and boundary counts
2. coordinates within a declared tolerance
3. cell types and ordered connectivity
4. face node order and owner/neighbour relationships
5. boundary patch/type assignments
6. interface and partition maps when exported
7. field association, component count, and values when fields are included

If an adapter cannot preserve a property, its capability declaration and test should state that limitation.

## 5. Suggested package layout

A diagnostic package could initially contain:

- `manifest.json`: schema version, dimension, counts, ID conventions, payload list, and capability flags
- `geometry.json` or a compact binary payload: coordinates and coordinate metadata
- `cells.json` or a compact binary payload: cell types and connectivity
- `faces.json` or a compact binary payload: ordered face connectivity and owner/neighbour
- `boundaries.json`: boundary regions, types, and face membership
- `interfaces.json`: optional interface and partition maps
- `fields/`: optional field payloads and association metadata

JSON is useful for inspecting small examples and debugging. For large production meshes, use a binary payload with explicit integer widths, endianness, array lengths, and schema version. The manifest can remain human-readable. Avoid inventing a large custom format before the semantic model and round-trip tests are stable.

## 6. Visualization choices

Expose visualization as an adapter choice, not as a responsibility of the mesh core:

- **Text summary:** fast diagnostics and CI artifacts.
- **VTK/ParaView:** general-purpose interactive 2D/3D visualization.
- **Python/PyVista:** scripted inspection, regression images, and notebooks.
- **Python/Matplotlib:** lightweight 1D curves and 2D structured/simple unstructured cases.
- **CGNS/Tecplot workflows:** interoperability with existing CFD tools.

A visualization export is not automatically a lossless interchange export. Keep those contracts separate.

## 7. Incremental implementation plan

1. Define a read-only mesh summary view and test its counts against a tiny hand-built mesh.
2. Define an L1 snapshot with explicit ID conventions and validate all connectivity ranges.
3. Add L1 serialization and a read-back test.
4. Add L2 face and boundary semantics; test orientation, owner/neighbour, and boundary membership.
5. Add one visualization adapter, preferably VTK XML, and test that it opens or parses with expected point/cell counts.
6. Add L3 partition/interface maps only after the non-partitioned round trip is stable.
7. Add optional fields and derived-data rebuilding.
8. Add capability reports for each adapter so unsupported information is visible.

Each step should be a small commit with a focused unit test. A test that is difficult to write without constructing the whole solver is evidence that the read-only view or snapshot boundary may be missing.

## 8. Current architecture constraints

The existing topology containers expose mutable arrays, and some derived adjacency is refreshed explicitly rather than automatically. Export must therefore either:

- run a validator and explicitly rebuild required derived data before taking a snapshot, or
- reject an inconsistent mesh with a useful diagnostic.

Do not claim snapshot consistency merely because an export function completed. A later refactor can centralize topology mutation and cache invalidation; that should be a separate, test-driven change.

## 9. Acceptance criteria

The first useful milestone is not a polished viewer. It is a tiny mesh that can be:

1. created with known geometry and topology,
2. exported to the diagnostic package,
3. read into an independent snapshot,
4. checked against the original semantic invariants,
5. exported to a visualization format without changing the canonical snapshot.

This milestone establishes a stable foundation for multiple viewers and interchange formats.
