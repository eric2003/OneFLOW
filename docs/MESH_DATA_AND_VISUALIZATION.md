# Mesh Data, Reconstruction, and Visualization Architecture

## Status

This document defines an incremental target architecture. It does not replace current mesh classes or file formats in one step. Each implementation step must preserve existing behavior and add focused tests.

## 1. Design goals

OneFLOW should be able to:
- inspect a mesh without loading every derived structure;
- export a useful summary and progressively richer mesh information;
- reconstruct a mesh from an exported data package;
- convert between the in-memory model and external formats;
- visualize the same mesh through external tools or a small Python viewer;
- test geometry and topology contracts independently of a full CFD solver.

The in-memory model, interchange schema, file writers, and visualization backends are separate responsibilities. A visualization format must not become the authoritative mesh model.

## 2. Information levels

The levels are cumulative in meaning, but a package may include only the levels it needs.

| Level | Data | Typical use | Reconstruction guarantee |
| --- | --- | --- | --- |
| L0 Summary | Dimension, counts, element-type counts, bounds, region count, schema version | Fast inspection and diagnostics | No mesh reconstruction |
| L1 Geometry | Node coordinates, coordinate dimension, coordinate system and units when known | Shape display and geometric queries | Geometry only |
| L2 Topology | Element types and ordered element-to-node connectivity; canonical face connectivity and owner/neighbour data when required by the selected mesh model | Mesh skeleton, adjacency, topology checks | Geometry and topology |
| L3 Semantics | Zone and region IDs, boundary patches, boundary names/types, physical tags | CFD setup and boundary visualization | Solver-relevant semantics, if all required fields are present |
| L4 Partition | Partition ownership, global/local ID maps, interface face/zone maps, ghost/halo metadata when required | Distributed mesh reconstruction and interface checks | Partition structure, if maps are complete |
| L5 Provenance | Generator inputs, layout parameters, conversion options, source-file fingerprints, software/schema versions | Reproducibility and mesh comparison | Supports regeneration; does not guarantee bitwise-identical output across versions |

L0 is derived from other levels and should not be treated as the canonical mesh. L1 and L2 are the minimum useful interchange package for general mesh reconstruction. L3 and L4 are optional only when the consumer does not need their semantics.

A package must declare which levels are present and which guarantees it provides. Missing data must not be silently guessed.

## 3. Canonical versus derived data

Prefer one authoritative representation for each fact.

- Node coordinates and connectivity are canonical mesh data.
- For a face-based finite-volume representation, ordered face-node connectivity and owner/neighbour relationships preserve orientation and flux semantics.
- Boundary patches and physical labels are canonical semantic data, not display colors.
- Cell-to-cell adjacency, search structures, centers, normals, areas, volumes, and other metrics should be derived or cached with explicit invalidation rules unless a format contract requires them as source data.
- Ghost entities and partition interfaces belong to the distribution/execution layer. They must not be confused with physical cells in the source mesh.
- Generator parameters describe provenance; they are not a substitute for the generated mesh.

The exact canonical topology should remain compatible with OneFLOW's actual mesh model and supported importers. Do not introduce duplicate mutable connectivity merely to make a writer convenient.

## 4. Interchange package

The first implementation should use a documented, versioned package rather than bind the design to one visualization format.

A package can contain:

- manifest.json: schema version, available levels, dimension, counts, coordinate units, source format, and file checksums;
- geometry: node coordinates;
- topology: element types, connectivity, and required face relationships;
- semantics: zones, boundary patches, names, and physical tags;
- partition: partition IDs and global/local/interface mappings;
- provenance: source inputs and generation/conversion settings.

Large numeric arrays should be stored in a compact, typed representation rather than expanded into verbose JSON. JSON is suitable for the manifest and small metadata. Binary or text array encodings should be chosen behind a format adapter. All identifiers must state whether they are zero-based or one-based; the schema must not rely on undocumented conventions.

A reader should validate schema version, array lengths, identifier ranges, element types, and cross-references before constructing live mesh objects. Failure must not leave a partially reconstructed mesh visible to callers.

## 5. Reconstruction and validation pipeline

Use a pipeline with explicit boundaries:

1. MeshInspector derives L0 counts and bounds from a read-only mesh view.
2. MeshSnapshot captures the selected canonical geometry, topology, and semantic data without owning solver state or transient caches.
3. MeshPackageWriter writes a versioned package through a format-specific adapter.
4. MeshPackageReader parses data into a neutral representation and validates it before committing a result.
5. MeshReconstructor builds the requested in-memory mesh representation.
6. MeshValidator checks topology, geometry, semantics, and partition invariants.
7. A visualization adapter consumes a read-only snapshot and does not mutate the mesh.

These are responsibility names, not a requirement to create all classes immediately. Begin with small pure functions or value types if they are easier to test. Avoid a universal manager that owns parsing, mesh mutation, validation, and rendering.

## 6. Visualization choices

The same canonical snapshot should support multiple outputs.

| Target | Best use | Notes |
| --- | --- | --- |
| VTK / VTU | General scientific visualization, mixed cell types, fields, and integration with ParaView | Strong default for visualization export |
| CGNS | CFD interchange, zones, boundary conditions, and flow solutions | Preserve CGNS-specific semantics in its adapter |
| SU2 mesh | CFD workflows that consume SU2 mesh input | Not a complete replacement for all OneFLOW metadata |
| Plot3D | Existing workflows and compatibility with established CFD tools | Keep legacy layout and format conventions explicit |
| Python with NumPy and Matplotlib | Lightweight 1D/2D plots, debugging, teaching, automated snapshots | Good initial in-repository viewer; 3D support is limited |
| Python with PyVista | Interactive 3D inspection, cell/point coloring, clipping, and screenshots | Optional dependency; can consume VTK-compatible data |
| Text / CSV summaries | CI diagnostics, small meshes, topology debugging | Useful for inspection, not a scalable full-mesh format |

Visualization modes should be selectable independently: points, surface edges, cell regions, boundary patches, partition IDs, face orientation, and scalar fields. A viewer should clearly distinguish physical boundaries from partition interfaces and ghost entities.

## 7. Test strategy

Tests should prove data meaning, not only that a file exists.

- Summary tests: counts, type histogram, and coordinate bounds for a tiny known mesh.
- Snapshot tests: the exported snapshot contains the expected ordered connectivity and labels.
- Reconstruction tests: export then import preserves counts, coordinates, connectivity, owner/neighbour pairs, boundary labels, and partition maps included in the selected levels.
- Negative tests: invalid IDs, mismatched array lengths, unsupported schema versions, and invalid cross-references are rejected without publishing a partial mesh.
- Geometry tests: positive cell measures and valid face orientation for supported examples.
- Partition tests: each local interface mapping resolves to the intended zone and face.
- Visualization tests: verify adapter output structure and metadata; keep image/pixel comparisons optional because rendering can vary by platform.

Round-trip equivalence must be defined semantically. Byte-for-byte equality is not required when writers legitimately reorder arrays, but identifiers and connectivity must be remapped consistently.

The small-mesh tests should not need a solver, a GUI, a global project, or external visualization software. If they do, treat that coupling as an architecture signal and reduce it where practical.

## 8. Incremental implementation plan

1. Keep the current topology model and add a read-only mesh summary API with unit tests.
2. Define a small neutral snapshot type for geometry and connectivity, with explicit index conventions.
3. Add validation for the snapshot independently of file I/O.
4. Implement a first package writer and reader, then round-trip tests.
5. Add L3 boundary/region metadata and test semantic preservation.
6. Add L4 partition mappings only after existing partition invariants are captured in tests.
7. Add VTK/VTU export and a minimal Python plotting script as consumers of the same snapshot.
8. Add further adapters only when a concrete workflow needs them.

At each step, inspect actual call sites and preserve existing CFD regressions. Do not claim a local test baseline still passes until the user or CI runs the tests after the change.
