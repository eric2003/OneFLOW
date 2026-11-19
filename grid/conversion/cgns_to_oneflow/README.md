# CGNS to OneFLOW conversion stage

This reusable Case is the CGNS input stage of the Mesh-to-CFD Validation
workflow. It consumes the CGNS artifact produced by the preceding
`grid/conversion/gmsh_to_cgns` Case and writes a native OneFLOW unstructured grid.

The workflow graph supplies the producer Case directory and relative CGNS
artifact path through `sourceGridCaseDir` and `sourceGridFileName`. It also
binds the output names for the native grid, boundary-region map, and grid
information report to the selected mesh Case. The checked-in parameters are
defaults only; the staged Case receives runtime artifact bindings and is
executed outside the source repository.

This is a conversion-stage template, not an STL-to-Gmsh geometry example.
It is deliberately not registered as an independent mesh-generation Case:
without an upstream CGNS artifact it has no standalone input.

The conversion preserves physical patch names and emits the boundary-region
map separately from the native grid. Downstream CFD consumes both artifacts.

When executed as part of Mesh-to-CFD Validation, the CGNS input is accompanied
by a checksum-verified ``.geometry-lineage.json`` record. After conversion, the
stage writes the inherited record beside the OneFLOW native grid and links it
to the exact CGNS input checksum. These sidecars are workflow provenance
metadata, not part of the native grid binary format.
