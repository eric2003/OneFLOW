# STL format inspection case and reusable readers

This directory has two related responsibilities: reusable STL readers/topology checks and a small standalone **STL inspection case**. Run this case first when bringing a new STL file into the Grid workflow; it checks that the file can be read and emits a compact report before any conversion or meshing is attempted.

## First-stage inspection case

From the repository root:

```sh
python grid/stl/script/inspect.py
```

To inspect another STL or choose an output location:

```sh
python grid/stl/script/inspect.py --input path/to/model.stl --output path/to/stl_report.json
```

The default input is `input/body.stl`; the default report is `results/stl_report.json`. Results are generated and ignored by Git. The report includes detected encoding (ASCII/Binary), triangle count, coordinate bounds, planarity, and—when the surface is a supported single planar component—boundary vertex/edge counts. Invalid or unsupported geometry returns a non-zero exit code with a diagnostic.

## Supported input

- ASCII STL: facet vertex records.
- Binary STL: 80-byte header, triangle count, and 50-byte triangle records.
- Finite coordinates are required. The first 2D conversion path accepts one connected, closed triangulated loop in a constant-Z XY plane.

## Current topology checks

The planar boundary pipeline welds coincident vertices using a caller-selected tolerance, counts triangle-edge incidence, extracts singly incident boundary edges, and rejects degenerate triangles, non-manifold edges, open/branched boundaries, multiple loops and self-intersecting boundaries.

STL defines geometry and triangle topology, not CFD boundary conditions. Boundary names therefore come from an explicit JSON sidecar in the consuming case. Physical group names are labels only; mapping them to CGNS or OneFLOW solver boundary types is a separate conversion responsibility.

## Tests

Run from the repository root:

```sh
python -m unittest grid.stl.test_stl_to_gmsh
```

The binary STL test builds its fixture in a temporary directory. Tests that generate reports or meshes use temporary directories and do not write generated files into this source directory.

## Scope

Arbitrarily oriented plane projection, multiple disconnected bodies, internal holes, non-planar surface repair, and a general-purpose production Mesh Core STL reader are not claimed by this first milestone.
