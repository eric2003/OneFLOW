# STL to Gmsh: native OneFLOW Cases

This directory contains native C++ STL-to-Gmsh Cases and independent Python
reference Cases. The native Cases run through the standard OneFLOW Grid task
and the C++ Geometry conversion pipeline.

| Implementation | Case | Layout | Nodes | Cells | Boundary elements |
|---|---|---:|---:|---:|---:|
| Native C++ | square/cpp | 1 x 1 | 4 | 1 QUAD_4 | 4 |
| Native C++ | square_2x2/cpp | 2 x 2 | 9 | 4 QUAD_4 | 8 |
| Native C++ | two_square_cells/cpp | 2 x 1 | 6 | 2 QUAD_4 | 6 |
| Native C++ | circle_external_flow/cpp | Unstructured | Gmsh-generated | TRI_3 | Boundary groups |\n| Native C++ | cube_single_hex/cpp | 1 x 1 x 1 | 8 | 1 HEXA_8 | 6 QUAD_4 |
| Python | square/python | 1 x 1 | 4 | 1 QUAD_4 | 4 |
| Python | two_square_cells/python | 2 x 1 | 6 | 2 QUAD_4 | 6 |
| Python | circle_external_flow/python | Unstructured | Gmsh-generated | TRI_3 | Boundary groups |
| Python | channel_internal_flow/python | 20 x 5 | 126 | 100 QUAD_4 | 50 |\n| Python | cube_single_hex/python | 1 x 1 x 1 | 8 | 1 HEXA_8 | 6 QUAD_4 |

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

The circle external-flow Cases use Inlet (left), Outlet (right), Farfield
(top and bottom, combined into one Physical Group), and Wall (the STL body).
The channel Case uses Inlet, Outlet, LowerWall and UpperWall. These names
identify patches; the Gmsh-to-CGNS adapter maps common aliases to CGNS
`BCType_t` values, and the CGNS-to-OneFLOW reader maps those broad categories
to OneFLOW boundary IDs. The original patch names remain intact. This is a
type-level mapping only: inlet/outlet state values and solver-specific wall
treatment require separate configuration.

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

## Shared output contract

Each implementation writes its mesh using the Case-specific canonical name: `square.msh`, `two_square_cells.msh`, `circle_external_flow.msh`, or `channel_internal_flow.msh`. Python declares this name in `expected.json` under `artifact.mesh_file`; its geometry and validation report use the same basename. In particular, the square outputs are `results/square.msh` and `results/square.stl-to-gmsh-report.json` for both C++ and Python. The producer language does not change the artifact naming contract, so downstream conversion and validation stages can consume either implementation without filename-specific branching.

## Python reference

Run the Python reference Cases through the staging orchestrator below. Do not
pass checked-in source Case directories to the Python runner directly, because
that would place generated files in the repository.

The Python Cases use the independent grid/stl/stl_to_gmsh.py implementation. They emit the same Case-specific mesh and report filenames as their native counterparts. The native C++ and Python `circle_external_flow` Cases use byte-identical STL input and the same domain, mesh-size and Physical Group settings. The native Case exercises OneFLOW's C++ STL reader, outer-domain/hole construction, Gmsh invocation and MSH validation; the Python Case is an independent implementation reference. The `circle_external_flow` Case treats a polygonal STL obstacle as an internal Wall loop and subtracts it from a rectangular outer domain. The `channel_internal_flow` Case treats the STL loop itself as the fluid domain, with inlet, outlet and wall groups; a transfinite 20 x 5 layout produces a reproducible quadrilateral channel mesh. The external case validates minimum mesh counts rather than exact counts because unstructured meshing can vary slightly between Gmsh versions. Both implementations also validate the generated coordinate bounds against the configured planar domain, reject degenerate triangles, and enforce the shared minimum normalized triangle quality (0.01) and minimum angle (1 degree). Python reports the measured minimum area, normalized quality, and minimum angle; native C++ reports the same quality thresholds alongside its MSH topology report. The two triangulations are not expected to have identical element counts or connectivity.

## Run all ten development Cases

After building OneFLOW and installing Gmsh, run the orchestrator from the
repository root:

    python grid/stl_to_gmsh/develop_test/run_cases.py

If the executable is not in the standard build/bin locations:

    python grid/stl_to_gmsh/develop_test/run_cases.py --oneflow-executable build/bin/Debug/OneFLOW.exe --gmsh gmsh

The orchestrator stages all ten Cases under
%USERPROFILE%/.oneflow/mesh-layout-test/grid/stl_to_gmsh on Windows or
$HOME/.oneflow/mesh-layout-test/grid/stl_to_gmsh on Linux. Staging paths mirror
the Case layout, for example square/cpp and square/python; each has its own
results directory. It recreates only the ten named Case staging directories
and leaves generated meshes, reports, and logs outside the repository source tree.

## STL-to-CFD pipeline

The optional dependency workflow uses the same square geometry and the
2 x 2 (four-cell) layout. It executes four standard OneFLOW Cases in order:

1. STL -> Gmsh: native C++ Grid task creates and validates results/square_2x2.msh.
2. Gmsh -> CGNS: the existing conversion Case consumes that MSH artifact and writes results/square.cgns.
3. Gmsh -> OneFLOW: the pipeline Case consumes the CGNS interchange artifact and writes
   Case-specific artifacts under ``results/``: ``<mesh-stem>.oneflow.grid``,
   ``<mesh-stem>_bcRegionMap.txt``, and ``<mesh-stem>.info.json``. For example,
   the cube prototype writes ``cube_single_hex_bcRegionMap.txt`` and
   ``cube_single_hex.info.json``. The boundary-map and grid-info filenames are
   supplied as runtime output bindings, so the conversion parameter file does
   not leak the square prototype names into other Cases. The CFD solver reads
   the boundary-map sidecar from the producing stage through the artifact binding.
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

The Python development unit tests live under `tests/python/grid`; the executable Cases themselves remain under `grid/stl_to_gmsh`.\n\nThe manually dispatched GitHub Actions workflow named OneFLOW External Tool
Integration runs the native C++ Gmsh integration tests, all ten development
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

    python -m unittest tests.python.grid.test_stl_to_gmsh_cases
    python tools/mesh_cases.py validate stl_to_gmsh_square
    python tools/mesh_cases.py validate stl_to_gmsh_square_2x2
    python tools/mesh_cases.py validate stl_to_gmsh_two_square_cells
    python tools/mesh_cases.py validate circle_external_flow/cpp

When Gmsh is available, the Python test suite executes the circle reference
Case and cube_single_hex in temporary directories, validating their MSH 4.1
element topology and Physical Groups. The native C++ circle integration test is enabled
by `ONEFLOW_ENABLE_GMSH_INTEGRATION_TESTS` and is exercised by the
`OneFLOW External Tool Integration` workflow.

Tests use temporary output directories. The structured native Cases support
one connected, axis-aligned rectangular XY boundary and QUAD_4 layouts. The
native circle external-flow Case additionally supports one closed planar STL
body loop inside a rectangular outer domain, producing a triangular fluid
region with a hole. The cube_single_hex pair is the only supported 3-D volume
prototype and is restricted to an axis-aligned cube. Arbitrary multiple body
loops, self-intersections, curved CAD entities, and general STL-to-volume
meshing are not claimed.

## Retired legacy entry point

The former root-level `script/run.py` workflow was removed because it wrote generated files into the source Case tree. Use the registered `circle_external_flow/cpp` and `circle_external_flow/python` Cases through the isolated development runner instead.

## 3D cube single-hexahedron prototype

The `cube_single_hex` pair is a deliberately restricted 3D STL-to-Gmsh prototype. Both implementations independently validate an axis-aligned, closed cube STL with equal edge lengths and exactly 12 triangles, eight bounding-box corners, two triangles per face, and a closed two-manifold edge incidence. They reconstruct the six planar faces and one volume, then ask Gmsh for a transfinite/recombined volume mesh.

The acceptance contract is exactly 8 nodes, 6 QUAD_4 boundary elements and 1 HEXA_8 volume cell. Physical Groups map XMin to `inflow`, XMax to `outflow`, YMin to `wall`, YMax to `farfield`, and both ZMin and ZMax to `symmetry`; `fluid_domain` is the volume group. These labels are subsequently mapped to CGNS BCType_t by the existing conversion layer.

This is not a general STL-to-hexahedral mesher. It accepts only the axis-aligned cube topology described above; arbitrary closed STL geometry requires a different volume-meshing strategy. Mesh-to-CFD Validation uses this case to verify STL/Gmsh/CGNS/OneFLOW topology and boundary handoff, and deliberately does not run the flat-plate CFD template on this one-cell cube.
