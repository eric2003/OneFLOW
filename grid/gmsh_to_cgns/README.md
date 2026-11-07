# Gmsh to CGNS bootstrap Case

This is an experimental, GUI-only manual Case for the Gmsh input adapter.
The primary sample is ASCII MSH 4.1; the C++ reader also retains the legacy
ASCII MSH 2.x path. Binary MSH and MSH 4.0 are not currently supported.

MSH 4.1 ASCII conversion writes a real CGNS unstructured Zone for supported
first-order top-dimensional cells. It also maps first-order elements from
the immediately lower dimension through Gmsh Entities and Physical Names
into CGNS boundary element sections and BC_t point lists. The conversion recognizes
explicit Physical Group names `inflow`, `outflow`, `symmetry`, and `farfield`, mapping them to the
corresponding CGNS `BCType_t` values. Unrecognized names remain
`BCTypeUserDefined`; names are never guessed from geometry or position. The writer reuses
OneFLOW's shared UnstructuredMeshData and UnstructuredCgnsZoneWriter, followed by the CgnsFactory / CgnsZbase file writer.
The bundled square sample writes four nodes, two TRI_3 cells, four BAR_2
boundary elements, and four named BC_t records.

## Parameters

The Case parameters are in `script/grid.txt`:

```cpp
string sourceGridType = "gmsh";
string targetGridType = "cgns";
```

## Current reader contract

- MSH 4.1 ASCII: entity-block node and element sections; sparse node tags;
  parametric node coordinates; common element types 1-31.
- MSH 2.x ASCII: legacy line-oriented nodes and elements.
- Reject binary input, unsupported MSH revisions, unsupported element types,
  duplicate/non-positive tags, and element references to unknown nodes.
- Node coordinates, ordered connectivity, sparse source tags, entity ownership,
  `$Entities`, and quoted `$PhysicalNames` are parsed into the GmshMesh value model.
- Entity-to-physical-group membership is retained and used for lower-dimensional
  boundary elements. Periodicity, partition and ghost metadata are not yet mapped.

## Manual execution

Run this Case only from the OneFLOW Tests GUI using the experimental
`Gmsh-to-CGNS Case (manual)` option. The GUI runner copies the Case to a
system temporary directory before invoking the standard Case entry point, so
the checked-in Case is not modified.

The current milestone is a first-order unstructured conversion path built
around a shared staging model and CGNS writer. The Gmsh reader does not pass
through SU2: it maps its GmshMesh value data into UnstructuredMeshData and
calls UnstructuredCgnsZoneWriter. The next engineering steps are to evolve
this limited staging model toward the fuller Mesh Core value contract, extend
three-dimensional and higher-order conversion coverage, and validate
round-trip semantics against additional real-world meshes. No Gmsh API is
used; the reader depends only on the MSH file format.


## MSH 4.1 reader progress

The C++ reader materializes node coordinates, source node and element IDs,
geometric entity ownership, physical names, entity-to-physical-group membership,
element types, and ordered connectivity. Quoted physical names may contain
spaces. The conversion path maps top-dimensional cells and first-order
codimension-one physical-group elements into UnstructuredMeshData and writes
the CGNS Unstructured Zone through UnstructuredCgnsZoneWriter. SU2 uses the
same writer through its own adapter; Gmsh does not use the SU2 parser or
Su2Grid as an intermediate mesh container. MSH 2.x remains on the validation-only
path; high-order cells and periodic/partition metadata are not yet converted.


## Read diagnostics

For MSH 4.1 ASCII input, the conversion entry point writes two companion
reports under the Case `results/` directory: `*.gmsh-read-report.json` for
machine-readable diagnostics and `*.gmsh-report.md` for a human-readable
English inspection report. Both reports are written to the Case-local `results/` directory, independently of the CGNS target filename. The Markdown report includes the overall summary,
physical groups, geometric entities, element-type counts, node and element
previews, coordinate bounds, TRI_3 area and XY orientation statistics, CGNS
output status, and current limitations. Markdown is the reviewable source
artifact; PDF can be generated from it when an archival or publication copy
is needed.

`READ_SUCCESS` means the reader accepted the file structure, preserved the
source node tags and ordered connectivity, and passed its structural checks.
It does not by itself guarantee that conversion succeeded. The companion
CGNS conversion test independently reopens the output and checks zone,
coordinate, cell-section, boundary-section and BC_t data.


## Current CGNS conversion scope

- Writes a standards-based unstructured CGNS Base and Zone for MSH 4.1 ASCII.
- Maps first-order LINE_2, TRI_3, QUAD_4, TETRA_4, HEXA_8, PENTA_6 and PYRA_5 top-dimensional cells.
- Preserves source node ordering and maps sparse node tags before CGNS one-based numbering.
- Writes coordinates, top-dimensional cell sections, supported codimension-one boundary sections, BC_t point lists, and a read-report sidecar.
- Maps the explicit names inflow, outflow, symmetry, and farfield to CGNS BCInflow, BCOutflow, BCSymmetryPlane, and BCFarfield. Other names remain BCTypeUserDefined.
- MSH 2.x is parsed but not yet materialized for CGNS output. High-order cells, periodicity, and partition metadata are not yet mapped.


## Boundary semantic mapping

The current two-dimensional square fixture defines four separate curve
Physical Groups. Their names are the explicit semantic contract consumed by
the conversion adapter:

| Physical Group | Geometric edge | CGNS BCType_t |
|---|---|---|
| `inflow` | Left | `BCInflow` |
| `outflow` | Right | `BCOutflow` |
| `symmetry` | Bottom | `BCSymmetryPlane` |
| `farfield` | Top | `BCFarfield` |

The converter does not infer these labels from coordinates. If a group has an
unrecognized name, its CGNS type remains `BCTypeUserDefined`. The original
Physical Group name is retained as the BC_t name. This is a format conversion
mapping only; flow values and solver-specific boundary state data are not
created by this mapping.


## Small-mesh example and test strategy

Keep this directory as the manual, end-to-end Case. Grow its examples as
small, focused MSH fixtures rather than large CFD meshes. A proposed layout is:

```text
grid/gmsh_to_cgns/
  grid/square.msh                 # current 2-D conversion + four BC semantics
  examples/
    read_only/                    # parser coverage; no CGNS writer required
    convert_2d/                   # TRI_3 and QUAD_4
    convert_3d/                   # TETRA_4, HEXA_8, PENTA_6, PYRA_5
    higher_order_read_only/       # TRI_6, TETRA_10, etc.; parse only initially
    invalid/                      # malformed references and unsupported cases
```

Use the smallest topology that proves each contract: one triangle for a basic
read, a four-node square for 2-D boundary semantics, one tetrahedron with four
triangular faces for 3-D boundary handling, and one cell per additional volume
family. High-order fixtures should initially verify parsing and ordered node
connectivity only; conversion tests must expect a clear unsupported-element
diagnostic until high-order CGNS mapping is implemented.

Separate the automated checks into three layers:

1. **Reader tests**: format/version, sparse node tags, entity ownership,
   Physical Group membership, element type and node count, ordered
   connectivity, and malformed-input diagnostics. These tests do not require
   CGNS.
2. **Conversion round-trip tests**: write a temporary CGNS file, reopen it
   through the CGNS API, and compare Base/Zone dimensions, coordinates,
   element sections, connectivity, BC_t names, BCType_t values, and point lists.
   Keep fixtures tiny and deterministic so these tests are suitable for CI.
3. **Manual Case tests**: exercise the complete Case entry point, result
   reports, and GUI workflow. Do not make large or experimental GUI Cases
   prerequisites for the fast parser unit tests.

Each fixture should state the exact capability it proves and its expected
node/element/section counts. Avoid one oversized mesh that mixes many
features: separate tiny fixtures make failures local and diagnostics clear.
## Interoperability route

The generic `GridConfig` conversion entry now supports both:

- `gmsh -> cgns`: write a CGNS interchange file.
- `gmsh -> oneflow`: compose the same CGNS writer with the CGNS reader and
  OneFLOW unstructured-grid construction. The intermediate file is retained
  as `results/<source-stem>.gmsh-to-oneflow.cgns`.

The composed route is intentionally not a second direct Gmsh runtime
converter. It exercises the CGNS serialization boundary and keeps one
authoritative Gmsh-to-CGNS mapping. It retains
`results/<source-stem>.gmsh-to-oneflow.cgns` as the interchange artifact and
writes the native OneFLOW grid to the configured `targetFile`; when no
native target is configured, the default is
`results/<source-stem>.oneflow.grid`. These outputs must remain separate:
the OneFLOW binary writer must never overwrite the CGNS input. The automated
test checks both files and independently inspects the OneFLOW unstructured
grid (nodes, cells, faces and coordinates). Its current input scope remains
MSH 4.1 ASCII and the first-order element families listed above. The
conversion matrix in `grid/interoperability/conversion_matrix.json` is the
central inventory of supported routes.
