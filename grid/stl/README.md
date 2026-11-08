# STL support in Grid

This directory contains reusable STL input and planar-topology utilities used by Grid conversion workflows. It is not a GUI feature and does not depend on the OneFLOW desktop application.

## Supported input

- ASCII STL: facet vertex records.
- Binary STL: 80-byte header, triangle count, and 50-byte triangle records.
- Finite coordinates are required. The current 2D conversion path accepts one connected, closed triangulated loop in a constant-Z XY plane.

## Current topology checks

The planar boundary pipeline welds coincident vertices using a caller-selected tolerance, counts triangle-edge incidence, extracts singly incident boundary edges, and rejects degenerate triangles, non-manifold edges, open/branched boundaries, multiple loops and self-intersecting boundaries.

STL defines geometry and triangle topology, not CFD boundary conditions. Boundary names therefore come from an explicit JSON sidecar in the consuming case. Physical group names are labels only; mapping them to CGNS or OneFLOW solver boundary types is a separate conversion responsibility.

## Tests

Run from the repository root:

```sh
python -m unittest grid.stl.test_stl_to_gmsh
```

The binary STL test builds its fixture in a temporary directory. No test writes generated files into this source directory.

## Scope

Arbitrarily oriented plane projection, multiple disconnected bodies, internal holes, non-planar surface repair, and a general-purpose production Mesh Core STL reader are not claimed by this first milestone.
