# STL to Gmsh: native OneFLOW Cases

This directory contains native C++ STL-to-Gmsh Cases and independent Python
reference Cases. The native Cases run through the standard OneFLOW Grid task
and the C++ Geometry conversion pipeline.

| Implementation | Case | Layout | Nodes | 2D cells | Boundary LINE_2 |
|---|---|---:|---:|---:|---:|
| Native C++ | square/cpp | 1 x 1 | 4 | 1 QUAD_4 | 4 |
| Native C++ | square_2x2/cpp | 2 x 2 | 9 | 4 QUAD_4 | 8 |
| Native C++ | two_square_cells/cpp | 2 x 1 | 6 | 2 QUAD_4 | 6 |
| Native C++ | circle_external_flow/cpp | Unstructured | Gmsh-generated | TRI_3 | Boundary groups |
| Python | square/python | 1 x 1 | 4 | 1 QUAD_4 | 4 |
| Python | two_square_cells/python | 2 x 1 | 6 | 2 QUAD_4 | 6 |
| Python | circle_external_flow/python | Unstructured | Gmsh-generated | TRI_3 | Boundary groups |
| Python | channel_internal_flow/python | 20 x 5 | 126 | 100 QUAD_4 | 50 |

For the structured square Cases, the STL describes the outer domain boundary.
For circle_external_flow, the STL describes the inner body boundary and the
rectangular outer domain is supplied by the layout. In every case, STL
triangle diagonals are triangulation edges, not CFD cell interfaces. The
square and square_2x2 Cases use byte-identical STL geometry: only
transfinite_cells changes from 1 x 1 to 2 x 2. The two_square_cells Case
deliberately uses a different, 2 x 1 rectangle; its split comes from
transfinite_cells.x = 2, not from STL triangles.

## Boundary contract

The square mesh Cases use four outer Physical Curves: inflow (left), outflow
(right), symmetry (bottom), and farfield (top). Their 2-D region is
fluid_domain.

The circle external-flow Cases use inflow (left), outflow (right), farfield
(top and bottom, combined into one Physical Group), and wall (the STL body).
Their fluid region is named `fluid_domain`. These are geometry labels;
Gmsh-to-CGNS maps recognized names to CGNS boundary types, while solver
boundary-state values remain a separate concern.

## Native C++ execution

Make Gmsh available on PATH, or change gmshExecutable in script/grid.txt.

    python grid/stl_to_gmsh/develop_test/run_cases.py --oneflow-executable build/bin/Debug/OneFLOW.exe --gmsh gmsh

The runner copies each Case to the external staging root before execution. Do
not execute a checked-in Case directly: relative output paths would create
runtime files in the source tree.

Linux uses the same Case paths with build/bin/OneFLOW. The C++ path reads STL,
recovers its planar boundary, writes .geo, invokes Gmsh, reads MSH 4.1 through
the existing Gmsh reader, and validates node/element counts and Physical
Group names. Outputs are written only under each Case results/ directory.

## Python reference

    python grid/stl_to_gmsh/run_python_case.py --case-directory grid/stl_to_gmsh/square/python
    python grid/stl_to_gmsh/run_python_case.py --case-directory grid/stl_to_gmsh/two_square_cells/python
    python grid/stl_to_gmsh/run_python_case.py --case-directory grid/stl_to_gmsh/circle_external_flow/python
    python grid/stl_to_gmsh/run_python_case.py --case-directory grid/stl_to_gmsh/channel_internal_flow/python

The Python Cases use the independent grid/stl/stl_to_gmsh.py implementation. The native C++ and Python `circle_external_flow` Cases use byte-identical STL input and the same domain, mesh-size and Physical Group settings. The native Case exercises OneFLOW's C++ STL reader, outer-domain/hole construction, Gmsh invocation and MSH validation; the Python Case is an independent implementation reference. The `circle_external_flow` Case treats a polygonal STL obstacle as an internal Wall loop and subtracts it from a rectangular outer domain. The `channel_internal_flow` Case treats the STL loop itself as the fluid domain, with inlet, outlet and wall groups; a transfinite 20 x 5 layout produces a reproducible quadrilateral channel mesh. The external case validates minimum mesh counts rather than exact counts because unstructured meshing can vary slightly between Gmsh versions. Both implementations also validate the generated coordinate bounds against the configured planar domain, reject degenerate triangles, and enforce the shared minimum normalized triangle quality (0.01) and minimum angle (1 degree). Python reports the measured minimum area, normalized quality, and minimum angle; native C++ reports the same quality thresholds alongside its MSH topology report. The two triangulations are not expected to have identical element counts or connectivity.

## Run all eight development Cases

After building OneFLOW and installing Gmsh, run the orchestrator from the
repository root:

    python grid/stl_to_gmsh/develop_test/run_cases.py

If the executable is not in the standard build/bin locations:

    python grid/stl_to_gmsh/develop_test/run_cases.py --oneflow-executable build/bin/Debug/OneFLOW.exe --gmsh gmsh

The orchestrator stages all eight Cases under
%USERPROFILE%/.oneflow/mesh-layout-test/grid/stl_to_gmsh on Windows or
$HOME/.oneflow/mesh-layout-test/grid/stl_to_gmsh on Linux. It recreates only
the eight named Case staging directories and leaves generated meshes, reports,
and logs outside the repository source tree.

## STL-to-CFD pipeline

The optional dependency workflow uses the same square geometry and the
2 x 2 (four-cell) layout. It executes four standard OneFLOW Cases in order:

1. STL -> Gmsh: native C++ Grid task creates and validates results/square_2x2.msh.
2. Gmsh -> CGNS: the existing conversion Case consumes that MSH artifact and writes results/square.cgns.
3. Gmsh -> OneFLOW: the pipeline Case consumes the CGNS interchange artifact and writes
   results/square.oneflow.grid, results/bcRegionMap.txt and
   results/square.info.json. The boundary map is the sidecar that the CFD
   solver reads for boundary-region setup; the info JSON records native-grid
   metadata for inspection.
4. CFD Solve: the dedicated cfd_uns2dslau2 pipeline Case consumes the native grid and must produce residual, restart-flow and visualization outputs.

CGNS remains an explicit, validated intermediate artifact; the end-to-end
consumer stage is named `gmsh_to_oneflow`. The workflow reuses the Workbench
dependency graph, artifact bindings and isolated case staging. It does not duplicate generated files into each
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

Native and Python circle implementations are staged together under
`circle_external_flow/cpp` and `circle_external_flow/python`, respectively.
Every process runs with its staged Case directory as the working directory;
no development Case is executed from the repository root.

The manually dispatched GitHub Actions workflow named OneFLOW External Tool
Integration runs the native C++ Gmsh integration tests, all eight development
Cases (Native C++ and Python Reference), and this dependency pipeline on
Windows. It uploads the isolated Case outputs and pipeline evidence as one
artifact.

Each Case declares its own output-validation rule in the workflow JSON, so
validation behavior is data-driven rather than selected by Case ID. The
runner snapshots the workflow JSON and writes one validation.json per
validated stage. The pipeline writes pipeline-report.json, environment.json,
run-manifest.json and workflow.json beside its staged Cases. Every executed
Case also receives results/runtime-environment.json beside its log. It independently checks the MSH
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
    python tools/mesh_cases.py validate circle_external_flow/cpp

When Gmsh is available, the Python test suite also executes the circle
reference Case in a temporary directory and validates its MSH 4.1 triangle
cells and Physical Groups. The native C++ circle integration test is enabled
by `ONEFLOW_ENABLE_GMSH_INTEGRATION_TESTS` and is exercised by the
`OneFLOW External Tool Integration` workflow.

Tests use temporary output directories. The structured native Cases support
one connected, axis-aligned rectangular XY boundary and QUAD_4 layouts. The
native circle external-flow Case additionally supports one closed planar STL
body loop inside a rectangular outer domain, producing a triangular fluid
region with a hole. Arbitrary multiple body loops, self-intersections,
curved CAD entities and 3-D volume meshing are not claimed.

## Legacy external-flow workflow

The original external-flow preprocessing entry point remains available for backward compatibility. The native `circle_external_flow` Case and its `circle_external_flow/python` reference are now the preferred executable examples for the external-flow domain model. It reads input/body.stl and applies
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
