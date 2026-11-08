# STL to Gmsh: native OneFLOW Cases

This directory contains native C++ STL-to-Gmsh Cases and independent Python
reference Cases. The native Cases run through the standard OneFLOW Grid task
and the C++ Geometry conversion pipeline.

| Implementation | Case | Transfinite cells | Nodes | QUAD_4 | Boundary LINE_2 |
|---|---|---:|---:|---:|---:|
| Native C++ | square | 1 x 1 | 4 | 1 | 4 |
| Native C++ | square_2x2 | 2 x 2 | 9 | 4 | 8 |
| Native C++ | two_square_cells | 2 x 1 | 6 | 2 | 6 |
| Python | python/square | 1 x 1 | 4 | 1 | 4 |
| Python | python/two_square_cells | 2 x 1 | 6 | 2 | 6 |

The STL files describe the outer boundary. Their triangle diagonals are
triangulation edges, not CFD cell interfaces. The square and square_2x2 Cases
use byte-identical STL geometry: only transfinite_cells changes from 1 x 1
to 2 x 2. The two_square_cells Case deliberately uses a different, 2 x 1
rectangle; its split comes from transfinite_cells.x = 2, not from STL triangles.

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

## Run all five development Cases

After building OneFLOW and installing Gmsh, run the orchestrator from the
repository root:

    python grid/stl_to_gmsh/develop_test/run_cases.py

If the executable is not in the standard build/bin locations:

    python grid/stl_to_gmsh/develop_test/run_cases.py --oneflow-executable build/bin/Debug/OneFLOW.exe --gmsh gmsh

The orchestrator stages all four Cases under
%USERPROFILE%/.oneflow/mesh-layout-test/grid/stl_to_gmsh on Windows or
$HOME/.oneflow/mesh-layout-test/grid/stl_to_gmsh on Linux. It recreates only
the five named Case staging directories and leaves generated meshes, reports,
and logs outside the repository source tree.

## STL-to-CFD pipeline

The optional dependency workflow uses the same square geometry and the
2 x 2 (four-cell) layout. It executes four standard OneFLOW Cases in order:

1. STL -> Gmsh: native C++ Grid task creates and validates results/square_2x2.msh.
2. Gmsh -> CGNS: the existing conversion Case consumes that MSH artifact and writes results/square.cgns.
3. CGNS -> OneFLOW: the pipeline CGNS consumer Case writes results/square.oneflow.grid.
4. CFD Solve: the dedicated cfd_uns2dslau2 pipeline Case consumes the native grid and must produce residual, restart-flow and visualization outputs.

The workflow reuses the Workbench dependency graph, artifact bindings and
isolated case staging. It does not duplicate generated files into each
consumer: sourceGridCaseDir plus the relative artifact path identifies the
producer output. Each staged Case and its logs live below the user's
.oneflow/mesh-layout-test/grid/stl_to_gmsh/pipelines directory. The source
Case tree and its parameter files are never modified.

Run only the end-to-end pipeline:

    python grid/stl_to_gmsh/develop_test/run_cases.py --pipeline-only

A different workflow definition can be supplied without changing the runner:

    python grid/stl_to_gmsh/develop_test/run_cases.py --pipeline-only --workflow path/to/workflow.json

Run the standalone geometry Cases first, followed by the pipeline:

    python grid/stl_to_gmsh/develop_test/run_cases.py --pipeline

Each Case declares its own output-validation rule in the workflow JSON, so
validation behavior is data-driven rather than selected by Case ID. The
runner snapshots the workflow JSON and writes one validation.json per
validated stage. The pipeline writes pipeline-report.json, environment.json
and run-manifest.json beside its staged Cases. It independently checks the MSH
4.1 node/line/quad counts and Physical Groups before allowing CGNS conversion
to proceed. It also requires non-empty residual, restart and visualization
outputs after CFD execution.

This first CFD stage is a runtime smoke test: it proves that the generated
mesh can be consumed and a solver can execute and write results. It is not
yet a claim of numerical accuracy, convergence, or grid-independence. Those
require separate residual/physical-state checks and analytical or manufactured
solution comparisons.

## Development tests

    python -m unittest grid.stl_to_gmsh.develop_test.test_cases
    python tools/mesh_cases.py validate stl_to_gmsh_square
    python tools/mesh_cases.py validate stl_to_gmsh_square_2x2
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
