# Gmsh to CGNS bootstrap Case

This is an experimental, GUI-only manual Case for the Gmsh input adapter.
The primary sample is ASCII MSH 4.1; the C++ reader also retains the legacy
ASCII MSH 2.x path. Binary MSH and MSH 4.0 are not currently supported.

The current output writer still emits a text placeholder with a `.cgns`
suffix. It is deliberately marked `NOT_A_CGNS_FILE`; CGNS tools must not
open or consume it. This Case validates input parsing and the normal
parameter-driven conversion entry point, not CGNS interoperability.

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
