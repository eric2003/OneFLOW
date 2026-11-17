# Conversion routes

A route documents a supported source-to-target workflow, its configuration,
input scope, intermediate artifacts and verification evidence.

Route-specific Cases belong here as they are migrated. Existing manual Cases
remain at their current paths until migration can be done without breaking
GUI and test references. Shared source fixtures belong under the format
owner, not inside each route.

The current Gmsh -> OneFLOW convenience route is implemented in
`GridConversion.cpp` and composes Gmsh -> CGNS -> OneFLOW. Its executable
round-trip test is `GmshCgnsConversionTest.ConvertsGmshThroughCgnsIntoOneFlowGrid`.
