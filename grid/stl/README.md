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

## Development test cases

STL reader compatibility cases are registered in `grid/case_catalog.json` under the `STL` group and appear in **Select Development Test Cases**. Each Python reference case lives under `grid/stl/python/<case>` and supplies `input/model.stl` and `expected.json`; the shared runner writes `results/stl_report.json` in the staged case workspace. Native OneFLOW cases live directly under `grid/stl/<case>` and run through the standard OneFLOW runtime. The native set currently includes `ascii_basic`, `ascii_crlf`, `ascii_scientific_notation`, `ascii_single_triangle`, `binary_basic`, `binary_header_solid`, `non_planar_surface`, and `weld_tolerance`.

The collection covers ASCII and binary STL, binary headers beginning with `solid`, CRLF text, scientific notation, native folded-surface conversion, configured vertex welding tolerance, degenerate triangles, truncated binary data, and Python planar-boundary diagnostics. Expected diagnostic failures count as passing tests when the observed diagnostic matches the case expectation.

Python development cases do not require a OneFLOW executable. In the Tests dialog, **Run Development Cases** enables native C++ cases and **Python reference cases** enables Python cases. Either option can run independently; enabling both runs both kinds in one staged batch. The selector shows only enabled kinds (and defaults to the native C++ catalog for browsing when neither run option is enabled). Selections are retained independently by execution kind. Python reference results do not substitute for native C++ evidence. A native Case is currently marked successful from its process exit status; its JSON mesh report is produced for inspection, but expected node/cell counts are not yet enforced by the GUI runner.

## STL to Gmsh generation

The first two-dimensional meshing workflow is maintained in
`grid/stl/stl_to_gmsh/`. It uses the shared `stl_to_gmsh.py` reader and
boundary reconstruction code to create Gmsh `.geo` geometry and MSH 4.1
ASCII meshes. It includes one unit-square QUAD_4 cell and two adjacent
unit-square QUAD_4 cells. The four outer Physical Curves are explicitly
named `inflow`, `outflow`, `symmetry`, and `farfield`; the fluid region
is `fluid_domain`. See that directory's README for staged execution,
expected node/cell counts, and the next Gmsh-to-CGNS handoff.

The cell layout is declared by `transfinite_cells.x` and
`transfinite_cells.y`, not inferred from STL triangle diagonals. Gmsh output
is MSH 4.1 ASCII so it can be consumed by the current C++ Gmsh reader.

## Tests

Run from the repository root:

```sh
python -m unittest grid.stl.test_stl_to_gmsh
python tools/test_mesh_cases.py validate
```

The binary STL test builds its fixture in a temporary directory. Tests that generate reports or meshes use temporary directories and do not write generated files into this source directory.

## Scope

Arbitrarily oriented plane projection, multiple disconnected bodies, internal holes, non-planar surface repair, and a general-purpose production Mesh Core STL reader are not claimed by this first milestone.


## Case directory convention

- Native OneFLOW case: `grid/stl/<case>`, catalog ID `stl_<case>`.
- Python reference case: `grid/stl/python/<case>`, catalog ID `stl_python_<case>`.

The case name describes the input variant (for example, `ascii_crlf`); the
directory and execution metadata identify the implementation path. The
catalog requires globally unique IDs, so Python references use the
`stl_python_` namespace.
