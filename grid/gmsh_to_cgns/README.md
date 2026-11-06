# Gmsh to CGNS bootstrap Case

This is an experimental, GUI-only manual Case for the Gmsh input adapter.
The primary sample is ASCII MSH 4.1; the C++ reader also retains the legacy
ASCII MSH 2.x path. Binary MSH and MSH 4.0 are not currently supported.

MSH 4.1 ASCII conversion now writes a real CGNS unstructured Zone for
first-order top-dimensional cells. The writer reuses OneFLOW's
Su2Grid::FillSU2CgnsZone() data-population path and CgnsFactory /
CgnsZbase serialization path. The bundled square sample writes two TRI_3
cells and four nodes. Lower-dimensional Gmsh elements are not yet written as
CGNS boundary sections; Physical Groups and Entities are the next mapping
stage.

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
- Node coordinates and element connectivity are parsed and validated, but
  are not yet transferred into OneFLOW's canonical mesh model.
- `$Entities`, `$PhysicalNames`, physical groups, periodicity, partition and
  ghost metadata are not yet preserved. Their mapping is a required next stage.

## Manual execution

Run this Case only from the OneFLOW Tests GUI using the experimental
`Gmsh-to-CGNS Case (manual)` option. The GUI runner copies the Case to a
system temporary directory before invoking the standard Case entry point, so
the checked-in Case is not modified.

The current milestone is an input-format adapter prototype. The next
engineering steps are to introduce a format-neutral mesh snapshot, preserve
entity and physical-group semantics as named boundaries/regions, and then
implement a real CGNS writer with semantic round-trip validation. No Gmsh API
is used; the reader depends only on the MSH file format.


## MSH 4.1 reader progress

The C++ reader now materializes node coordinates, source node and element IDs, geometric entity ownership, element types, and ordered connectivity in a GmshMesh value model. A focused GTest covers the square-mesh shape and malformed references. The conversion entry point invokes this reader and cross-checks its counts. The MSH 4.1 first-order top-dimensional cell data is now mapped through the existing Su2Grid-to-CgnsZone population path and written through CgnsFactory. The current stage omits lower-dimensional boundary elements and Physical Groups; MSH 2.x remains on the validation-only path.


## Read diagnostics

For MSH 4.1 ASCII input, the conversion entry point now writes a companion
`*.gmsh-read-report.json` next to the configured output. The console reports
the parsed node/element totals and TRI_3 geometry summary. The report includes
source tag ranges, coordinate bounds, element-type counts, TRI_3 area and XY
orientation statistics, and explicit validation states.

`READ_SUCCESS` means the reader accepted the file structure, preserved the
source node tags and ordered connectivity, and passed its structural checks.
It does not mean that Physical Groups have been mapped or that a CGNS file
has been produced. The current `.cgns` output remains a clearly marked
non-CGNS placeholder until the existing OneFLOW CGNS Zone path is connected.


## Current CGNS conversion scope

- Writes a standards-based unstructured CGNS Base and Zone for MSH 4.1 ASCII.
- Maps first-order LINE_2, TRI_3, QUAD_4, TETRA_4, HEXA_8, PENTA_6 and PYRA_5 top-dimensional cells.
- Preserves source node ordering and maps sparse node tags before CGNS one-based numbering.
- Writes coordinates, cell sections and a read-report sidecar.
- Lower-dimensional elements, PhysicalNames, Entities, and CGNS BC_t records are not yet emitted. No Physical Group name is inferred as a solver BC.
- MSH 2.x is parsed but not yet materialized for CGNS output. High-order cells are not yet mapped.
