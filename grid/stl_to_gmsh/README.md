# STL to Gmsh: native OneFLOW Cases

This directory contains native C++ STL-to-Gmsh Cases and independent Python
reference Cases. The native Cases run through the standard OneFLOW Grid task
and the C++ Geometry conversion pipeline.

| Implementation | Case | Transfinite cells | Nodes | QUAD_4 | Boundary LINE_2 |
|---|---|---:|---:|---:|---:|
| Native C++ | square | 1 x 1 | 4 | 1 | 4 |
| Native C++ | two_square_cells | 2 x 1 | 6 | 2 | 6 |
| Python | python/square | 1 x 1 | 4 | 1 | 4 |
| Python | python/two_square_cells | 2 x 1 | 6 | 2 | 6 |

The STL files describe the outer boundary. Their triangle diagonals are
triangulation edges, not CFD cell interfaces. The two-cell split at x = 1
comes from transfinite_cells.x = 2, not from STL triangles.

## Boundary contract

The four external Physical Curves are inflow (left), outflow (right),
symmetry (bottom), and farfield (top). The 2-D region is fluid_domain.
These are Gmsh physical labels, not solver boundary-state definitions.

## Native C++ execution

Make Gmsh available on PATH, or change gmshExecutable in script/grid.txt.

    build\bin\Debug\OneFLOW.exe r grid/stl_to_gmsh/square
    build\bin\Debug\OneFLOW.exe r grid/stl_to_gmsh/two_square_cells

Linux uses the same Case paths with build/bin/OneFLOW. The C++ path reads STL,
recovers its planar boundary, writes .geo, invokes Gmsh, reads MSH 4.1 through
the existing Gmsh reader, and validates node/element counts and Physical
Group names. Outputs are written only under each Case results/ directory.

## Python reference

    python grid/stl_to_gmsh/python/run_case.py --case-directory grid/stl_to_gmsh/python/square
    python grid/stl_to_gmsh/python/run_case.py --case-directory grid/stl_to_gmsh/python/two_square_cells

The Python Cases use the independent grid/stl/stl_to_gmsh.py implementation.

## Development tests

    python -m unittest grid.stl_to_gmsh.develop_test.test_cases
    python tools/mesh_cases.py validate stl_to_gmsh_square
    python tools/mesh_cases.py validate stl_to_gmsh_two_square_cells

Tests use temporary output directories. The first native stage supports one
connected, axis-aligned rectangular XY boundary and structured QUAD_4 layouts.
Multiple loops, holes, curved/rotated boundaries, unstructured size fields,
and 3-D volume meshing are not claimed.

## Existing external-flow workflow (preserved)

The original external-flow preprocessing case remains available alongside
the native square Cases. It reads input/body.stl and applies
config/external_flow.json. The STL body is an internal Wall boundary; the
rectangular fluid domain has Inlet, Outlet, Symmetry, and Farfield Physical
Curve groups, and the fluid region is named Fluid. These names are Gmsh
physical labels, not automatic OneFLOW solver boundary-condition mappings.

Run the original Python workflow from the repository root:

    python grid/stl_to_gmsh/script/run.py

To also invoke Gmsh and write the legacy MSH 2.2 output:

    python grid/stl_to_gmsh/script/run.py --gmsh gmsh

Its generated files remain under grid/stl_to_gmsh/results/:
body.geo, body.msh (when Gmsh is supplied), and mesh_report.json.
The checked-in body STL and external_flow.json are inputs and are not modified.
