# Grid Data Contracts and Visualization Architecture

## Goals

Keep mesh data independent from solver startup, global configuration, GUI code, and file-format libraries. Each layer must be testable from a small in-memory fixture.

The first deliverable is the L0 summary contract and its JSON writer. A summary is descriptive metadata; it is not sufficient to reconstruct a mesh.

## Proposed export layout

Use a directory per exported mesh so related files can be moved together:

```text
<mesh-name>.mesh/
  manifest.json
  summary.json
  geometry.json
  topology.json
  boundary.json
  views/
    mesh.vtu
    mesh.dat
```

All JSON documents must declare a schema identifier and version. IDs in exported files must have a documented convention (prefer zero-based IDs in interchange files, with explicit conversion at OneFLOW boundaries).

### L0: summary.json

Purpose: quick inspection, validation reports, UI previews, and batch inventory.

Fields:
- spatial dimension
- node, face, cell, and boundary-face counts
- optional later additions: grid name, coordinate units, cell-type histogram, bounding box, source format

This level does not contain coordinates or connectivity and cannot reconstruct the grid.

### L1: geometry.json

Purpose: portable point coordinates.

Fields:
- coordinate component names and units
- ordered node coordinates
- coordinate system metadata when known

The node ordering is part of the contract. Coordinates alone cannot reconstruct cell shapes or boundaries.

### L2: topology.json

Purpose: reconstruct cells and faces.

Fields:
- node IDs belonging to each face or cell
- face-to-cell adjacency where available
- cell type for mixed-element meshes
- explicit ID base and array offsets for compact encodings

Topology plus geometry is the minimum needed to reconstruct a geometric mesh. Face-to-cell adjacency can be derived for some cell-only formats, but should not be assumed when exact OneFLOW face ordering matters.

### L3: boundary.json

Purpose: preserve boundary semantics, not just shape.

Fields:
- boundary face IDs or ranges
- boundary-condition type/name
- patch or zone identity
- periodic/interface pairing metadata when applicable

A geometry/topology round trip is not a full CFD round trip unless boundary conditions and required interface metadata are preserved.

### manifest.json

Purpose: package-level index and compatibility contract.

Fields:
- manifest schema and version
- references to the summary, geometry, topology, and boundary documents
- coordinate units, ID base, source format, and export tool version
- optional checksums for detecting mismatched files

The loader must reject unsupported schema versions and mismatched referenced files with actionable errors.

## Output and visualization adapters

Keep a canonical in-memory mesh contract separate from output adapters. The intended flow is:

```text
Grid / UnsGrid
    -> summary + geometry + topology + boundary views
    -> canonical export model
    -> JSON package | VTK VTU | CGNS | Tecplot | Plot3D
                         |
                         -> Python (PyVista / meshio / VTK)
                         -> external visualization tools
                         -> optional OneFLOW viewer
```

Recommended order:
1. Implement and test the L0 summary and JSON output.
2. Expose summary data from existing grid objects without copying or owning mesh storage.
3. Define independent geometry and topology views and validate them with tiny fixtures.
4. Add boundary/interface metadata and round-trip tests.
5. Add one mature visualization adapter first (VTK VTU is a practical choice), then add other formats behind the same adapter boundary.
6. Provide a small Python example that reads the exported interchange format and displays it with PyVista. Python visualization should consume the same exported data as external tools, not a private duplicate format.

Do not make the core mesh model depend on PyVista, VTK, JSON libraries, or a GUI. Those belong in adapters.

## Test architecture requirements

- L0 tests must link only the summary contract and a test framework.
- Geometry tests must construct coordinates directly without launching a solver or reading global configuration.
- Topology tests must use small hand-built meshes and verify node references, cell references, face ownership, and boundary/interior distinctions.
- Export tests must verify schema/version, deterministic ordering, ID-base conversion, and clear failure on invalid input.
- Round-trip tests must compare dimension, coordinates, cell types, connectivity, boundary labels, and interface relations, not only counts.
- At least one integration test should export a tiny mesh and read it with an independent consumer when practical.

A test that requires full application initialization to check a local topology invariant is evidence that the dependency boundary should be reconsidered.

## Current scope

The current `GridSummary` is only L0. `WriteGridSummaryJson` serializes that summary with schema `oneflow.grid-summary`, version 1. Neither currently reads from `Grid`/`UnsGrid`, serializes topology, or promises mesh reconstruction. Those are explicit follow-up steps rather than implicit guarantees.
