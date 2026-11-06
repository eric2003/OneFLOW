# Gmsh to CGNS bootstrap Case

This is an experimental, GUI-only manual Case for the Gmsh input adapter.
The primary sample is ASCII MSH 4.1; the C++ reader also retains the legacy
ASCII MSH 2.x path. Binary MSH and MSH 4.0 are not currently supported.

MSH 4.1 ASCII conversion writes a real CGNS unstructured Zone for supported
first-order top-dimensional cells. It also maps first-order elements from
the immediately lower dimension through Gmsh Entities and Physical Names
into CGNS boundary element sections and BC_t point lists. Physical group
names are retained with BCTypeUserDefined; they are not inferred to be
solver-specific wall, inlet, or outlet conditions. The writer reuses
OneFLOW's Su2Grid::FillSU2CgnsZone() and CgnsFactory / CgnsZbase paths.
The bundled square sample writes four nodes, two TRI_3 cells, four BAR_2
boundary elements, and one named boundary record.

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

The current milestone is an input-format adapter prototype. The next
engineering steps are to introduce a format-neutral mesh snapshot, preserve
entity and physical-group semantics as named boundaries/regions, and then
implement a real CGNS writer with semantic round-trip validation. No Gmsh API
is used; the reader depends only on the MSH file format.


## MSH 4.1 reader progress

The C++ reader materializes node coordinates, source node and element IDs,
geometric entity ownership, physical names, entity-to-physical-group membership,
element types, and ordered connectivity. Quoted physical names may contain
spaces. The conversion path maps top-dimensional cells and first-order
codimension-one physical-group elements through the existing
Su2Grid-to-CgnsZone population path. MSH 2.x remains on the validation-only
path; high-order cells and periodic/partition metadata are not yet converted.


## Read diagnostics

For MSH 4.1 ASCII input, the conversion entry point now writes a companion
`*.gmsh-read-report.json` next to the configured output. The console reports
the parsed node/element totals and TRI_3 geometry summary. The report includes
source tag ranges, coordinate bounds, element-type counts, TRI_3 area and XY
orientation statistics, and explicit validation states.

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
- Preserves physical group names as CGNS boundary labels with BCTypeUserDefined; no solver boundary condition is inferred.
- MSH 2.x is parsed but not yet materialized for CGNS output. High-order cells, periodicity, and partition metadata are not yet mapped.
