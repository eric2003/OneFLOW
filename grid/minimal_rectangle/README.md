# Minimal Rectangle Mesh Case

This case generates a small structured mesh with OneFLOW, exports inspectable
metadata and visualization files, and checks the result against explicit
geometric and topological expectations.

## Case parameters

The case definition is `case.json`. It records the generation test, layout,
coordinate representation, expected counts, and expected bounding box. It is
configuration, not a place where calculations run.

The geometry comes from
`tests/geometry/data/minimalRectangleLayout.txt`: a unit square with three
nodes along each edge. The expected mesh has 9 nodes, 4 quadrilateral cells,
12 unique edges, 8 boundary edges, and 4 internal edges.

## Run in an external workspace

Build the existing `GridLayoutGenerationTest` target first. Then run from the
repository root, supplying its executable and a workspace outside the checkout:

```powershell
.\grid\minimal_rectangle\run.ps1 -TestExecutable .\build\tests\geometry\GridLayoutGenerationTest.exe -WorkspaceDir D:\OneFLOWWorkspace\cases\minimal_rectangle
```

Adjust the executable path to match the local build. The runner sends the
generated OneFLOW grid and exported package to the external workspace; it does
not write generated artifacts into the repository case directory.

The runner checks summary counts and the three-axis bounding box. The package
contains `summary.json`, `geometry.json`, `topology.json`, `boundary.json`,
`manifest.json`, `views/mesh.svg`, and `views/mesh.vtu`. Open the SVG in a
browser or the VTU in ParaView.

## Scope

This is a concrete generation-and-inspection case, not a general mesh
framework. Keep its parameters here and its outputs in the external workspace.
As further cases are added, share validation or visualization code only where
their real workflows establish the same need.
