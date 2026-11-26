# STL to Gmsh: native OneFLOW Cases

## Direct 3-D STL volume mesh command

Generate a tetrahedral volume mesh from a closed STL surface:

```powershell
python grid/stl_to_gmsh/generate_volume_mesh.py --stl path/to/surface.stl --output "$HOME/.oneflow/mesh-generation/surface/volume.msh" --mesh-size 0.20
```

To allow a size interval, specify both bounds:

```powershell
python grid/stl_to_gmsh/generate_volume_mesh.py --stl path/to/surface.stl --output "$HOME/.oneflow/mesh-generation/surface/volume.msh" --mesh-size-min 0.10 --mesh-size-max 0.20
```

Use the same STL and vary `--mesh-size` (for example 0.40, 0.20, 0.10)
to compare coarse, medium and fine meshes. The command writes MSH, GEO and
JSON validation-report files. The terminal ends with a concise SUCCESS/FAILED
summary; the JSON report records stage-level PASS/FAIL/NOT_RUN/NOT_CHECKED
states and the process exit code. See [the run-result contract](../../docs/source/run_result_reporting.rst)
and
[`stl_to_gmsh_commands`](../../docs/source/stl_to_gmsh_commands.rst).


Generate parameterized geometry from the Mesh Generation GUI or CLI before
meshing:

.. code-block:: powershell

   python grid/stl_to_gmsh/generate_stl_geometry.py --geometry sphere --radius 1 --subdivisions 3 --output "$HOME/.oneflow/mesh-generation/generated_stl/sphere/sphere.stl"
   python grid/stl_to_gmsh/generate_stl_geometry.py --geometry prolate_spheroid --semi-axis-a 3 --semi-axis-b 0.5 --semi-axis-c 0.5 --output "$HOME/.oneflow/mesh-generation/generated_stl/prolate_spheroid/prolate_spheroid.stl"
   python grid/stl_to_gmsh/generate_stl_geometry.py --geometry cylinder_internal_flow --radius 0.5 --length 4 --circumferential-divisions 48 --axial-divisions 16 --output "$HOME/.oneflow/mesh-generation/generated_stl/cylinder_internal_flow/cylinder_internal_flow.stl"

Each command also writes a .geometry.json manifest. For
`cylinder_internal_flow`, the STL and checksum-bound manifest are imported as
one mesh project. The Gmsh adapter validates the closed STL topology and patch
indices, then creates an analytic OpenCASCADE cylinder with separate
`wall`, `inflow`, and `outflow` Physical Surfaces plus the
`fluid_domain` Physical Volume. This preserves the intended CFD boundary
semantics instead of treating the end caps as solid walls. See
docs/source/stl_geometry_generation.rst.

## Incremental Case lifecycle contracts

New geometry and mesh Cases are delivered through explicit lifecycle contracts,
not an assumed all-at-once pipeline. The contract records which geometry,
surface-region, volume-mesh, conversion, OneFLOW-ingestion, CFD-smoke, and
physical-validation stages are planned, implemented, and actually verified.
A future stage being planned does not fail the acceptance checks for the
current delivery scope.

The common contract and stage-state rules are documented in
[OneFLOW Case Lifecycle and Incremental Verification](../../docs/architecture/oneflow-case-lifecycle-contract.md);
the machine-readable schema is
[case-lifecycle-v1.schema.json](../../tools/schemas/case-lifecycle-v1.schema.json).

The 90-degree pipe bend is the first staged example:
`3d/pipe_bend/90_degree/automatic_tetrahedra/case-lifecycle.json`.
Its lifecycle contract records Python and native C++ geometry/region
implementations independently. The native C++ producer is executable through
`run_pipe_bend_cpp_case.py`; it compiles into the staged Case workspace,
writes the same STL, checksum-bound region manifest and geometry manifest,
then delegates volume meshing to the shared
`grid/stl_to_gmsh/generate_volume_mesh.py` service. For example:

```powershell
python grid/stl_to_gmsh/3d/pipe_bend/90_degree/automatic_tetrahedra/run_pipe_bend_cpp_case.py --case-directory "$HOME/.oneflow/grid-development/pipe-bend-cpp" --gmsh gmsh
```

This is a native C++ geometry implementation driven by a Python workflow
adapter; it is not a native OneFLOW solver Case. In the OneFLOW Tests dialog,
enable the Python / mesh producer Cases, select the pipe bend, and choose
`Native C++ geometry` in its Implementation column. The Python reference
remains the default. The selected implementation is remembered by the dialog
and both variants run through the same staged Case workspace and Gmsh service.
The current lifecycle verification scope remains the Python reference
milestone until the native runner and its generated MSH4 output have passed
the complete declared acceptance checks. Downstream CGNS conversion, OneFLOW
ingestion and CFD stages remain separate milestones.

This directory contains native C++ STL-to-Gmsh Cases and independent Python
reference Cases. The native Cases run through the standard OneFLOW Grid task
and the C++ Geometry conversion pipeline. The geometric definitions and STL
construction rules for the sphere, prolate spheroid and internal-flow
cylinder are documented in [STL geometry generation](../../docs/source/stl_geometry_generation.rst).

| Implementation | Case | Layout | Nodes | Cells | Boundary elements |
|---|---|---:|---:|---:|---:|
| Native C++ | 2d/square/single_quad/cpp | 1 x 1 | 4 | 1 QUAD_4 | 4 |
| Native C++ | 2d/square/two_triangles/cpp | 1 x 1 | 4 | 2 TRI_3 | 4 |
| Native C++ | 2d/square/automatic_triangles/cpp | Automatic | Gmsh-generated | TRI_3 | Subdivided boundary |
| Native C++ | 2d/square_2x2/cpp | 2 x 2 | 9 | 4 QUAD_4 | 8 |
| Native C++ | 2d/two_square_cells/cpp | 2 x 1 | 6 | 2 QUAD_4 | 6 |
| Native C++ | 2d/circle_external_flow/cpp | Unstructured | Gmsh-generated | TRI_3 | Boundary groups |
| Native C++ | 3d/cube/single_hex/cpp | 1 x 1 x 1 | 8 | 1 HEXA_8 | 6 QUAD_4 |
| Native C++ | 3d/cube/six_tetrahedra/cpp | Fixed | 8 | 6 TETRA_4 | 12 TRI_3 |
| Native C++ | 3d/cube/automatic_tetrahedra/cpp | Automatic | Gmsh-generated | TETRA_4 | Triangular boundary |
| Native C++ | 3d/sphere/automatic_tetrahedra/cpp | Automatic | Gmsh-generated | TETRA_4 | Closed wall surface |
| Native C++ | 3d/prolate_spheroid/automatic_tetrahedra/cpp | Automatic | Gmsh-generated | TETRA_4 | Closed wall surface |
| Native C++ | 3d/cylinder_internal_flow/automatic_tetrahedra/cpp | Automatic | Gmsh OCC cylinder | TETRA_4 | Wall, inflow and outflow |
| Python | 2d/square/single_quad/python | 1 x 1 | 4 | 1 QUAD_4 | 4 |
| Python | 2d/square/two_triangles/python | 1 x 1 | 4 | 2 TRI_3 | 4 |
| Python | 2d/square/automatic_triangles/python | Automatic | Gmsh-generated | TRI_3 | Subdivided boundary |
| Python | 2d/two_square_cells/python | 2 x 1 | 6 | 2 QUAD_4 | 6 |
| Python | 2d/circle_external_flow/python | Unstructured | Gmsh-generated | TRI_3 | Boundary groups |
| Python | 2d/channel_internal_flow/python | 20 x 5 | 126 | 100 QUAD_4 | 50 |
| Python | 3d/cube/single_hex/python | 1 x 1 x 1 | 8 | 1 HEXA_8 | 6 QUAD_4 |
| Python | 3d/cube/six_tetrahedra/python | Fixed | 8 | 6 TETRA_4 | 12 TRI_3 |
| Python | 3d/cube/automatic_tetrahedra/python | Automatic | Gmsh-generated | TETRA_4 | Triangular boundary |
| Python | 3d/sphere/automatic_tetrahedra/python | Automatic | Gmsh-generated | TETRA_4 | Closed wall surface |
| Python | 3d/prolate_spheroid/automatic_tetrahedra/python | Automatic | Gmsh-generated | TETRA_4 | Closed wall surface |
| Python | 3d/cylinder_internal_flow/automatic_tetrahedra/python | Automatic | Gmsh-generated | TETRA_4 | Wall, inflow and outflow |

For the square Cases, the STL describes the outer domain boundary. The
single-quad and two-triangle variants use the same STL and 1 x 1 transfinite
boundary layout: the former recombines the surface into one QUAD_4, while the
latter leaves it un-recombined and produces two TRI_3 cells.
For the 2D circle external-flow Case, the STL describes the inner body boundary and the
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

The unit-sphere and 6:1 prolate-spheroid Cases also run the complete Mesh-to-CFD
smoke chain. Their closed volumes have a wall-only exterior (no inlet or
outlet); the workflow reuses the laminar-flat-plate CFD template for common
solver settings and overrides the staged dimension to 3. This deliberately
exercises the 3-D OneFLOW runtime after mesh handoff. It is a pipeline smoke
test, not a physically representative external-flow calculation or an
accuracy benchmark.

## Native C++ execution

Make Gmsh available on PATH, or change gmshExecutable in script/grid.txt.

    python grid/stl_to_gmsh/run_cases.py --oneflow-executable build/bin/Debug/OneFLOW.exe --gmsh gmsh

The runner copies each Case to the external staging root before execution. Do
not execute a checked-in Case directly: relative output paths would create
runtime files in the source tree.

Linux uses the same Case paths with build/bin/OneFLOW. The C++ path reads STL,
recovers its planar boundary, writes .geo, invokes Gmsh, reads MSH 4.1 through
the existing Gmsh reader, and validates node/element counts and Physical
Group names. Outputs are written only under each Case results/ directory.

## Shared output contract

Each implementation writes its mesh using the Case-specific canonical name: `square.msh`, `two_square_cells.msh`, `circle_external_flow.msh`, or `channel_internal_flow.msh`. Python declares this name in `expected.json` under `artifact.mesh_file`; its geometry and validation report use the same basename. In particular, the fixed square outputs are `results/square.msh`, while the automatic pair uses `results/square_automatic_triangles.msh`; each report shares its mesh basename. The producer language does not change the artifact naming contract, so downstream conversion and validation stages can consume either implementation without filename-specific branching.

## Python reference

Run the Python reference Cases through the staging orchestrator below. Do not
pass checked-in source Case directories to the Python runner directly, because
that would place generated files in the repository.

The Python Cases use `grid/stl_to_gmsh/converter.py`, built on the shared STL reader (`grid/stl_import/reader.py`) and surface-topology utilities (`grid/stl_import/surface.py`). They emit the same Case-specific mesh and report filenames as their native counterparts. The native C++ and Python `2d/square/automatic_triangles` Cases use byte-identical STL input and the same automatic layout. They omit transfinite constraints and rely on mesh_size for Gmsh refinement. The native C++ and Python `2d/circle_external_flow` Cases use byte-identical STL input and the same domain, mesh-size and Physical Group settings. The native Case exercises OneFLOW's C++ STL reader, outer-domain/hole construction, Gmsh invocation and MSH validation; the Python Case is an independent implementation reference. The `2d/circle_external_flow` Case treats a polygonal STL obstacle as an internal Wall loop and subtracts it from a rectangular outer domain. The `2d/channel_internal_flow` Case treats the STL loop itself as the fluid domain, with inlet, outlet and wall groups; a transfinite 20 x 5 layout produces a reproducible quadrilateral channel mesh. The external case validates minimum mesh counts rather than exact counts because unstructured meshing can vary slightly between Gmsh versions. Both implementations also validate the generated coordinate bounds against the configured planar domain, reject degenerate triangles, and enforce the shared minimum normalized triangle quality (0.01) and minimum angle (1 degree). Python reports the measured minimum area, normalized quality, and minimum angle; native C++ reports the same quality thresholds alongside its MSH topology report. The two triangulations are not expected to have identical element counts or connectivity.

## Run all registered STL-to-Gmsh Cases

After building OneFLOW and installing Gmsh, run the orchestrator from the
repository root:

    python grid/stl_to_gmsh/run_cases.py

If the executable is not in the standard build/bin locations:

    python grid/stl_to_gmsh/run_cases.py --oneflow-executable build/bin/Debug/OneFLOW.exe --gmsh gmsh

The standalone producer-check runner stages the 22 fixed Cases plus six sphere mesh-size variants under
%USERPROFILE%/.oneflow/grid-development/stl-to-gmsh on Windows or
$HOME/.oneflow/grid-development/stl-to-gmsh on Linux. Staging paths mirror
the Case layout, for example 2d/square/single_quad/cpp and 2d/square/single_quad/python; each has its own
results directory. It recreates only the 27 named Case staging directories
and leaves generated meshes, reports, and logs outside the repository source tree.

## Mesh-to-CFD validation

End-to-end validation is provided by the shared Mesh-to-CFD Validation
Workbench in `tools/oneflow_gui/mesh_validation.py`, not by a second
STL-specific pipeline runner. Select any supported native C++ or Python
STL-to-Gmsh Case; the Workbench builds the same artifact-bound chain:

1. STL -> Gmsh: execute the selected producer and inspect its MSH 4.1 output.
2. Gmsh -> CGNS: use the shared `grid/conversion/gmsh_to_cgns` conversion Case.
3. CGNS -> OneFLOW: use the shared `grid/conversion/cgns_to_oneflow` conversion-stage
   template, binding native-grid, boundary-region-map and grid-info outputs
   to the selected mesh stem.
4. CFD Solve: stage the existing `test/cfd/plateuns2dslau2` consumer and
   apply runtime overrides for a short smoke test. No STL-specific duplicate
   solver template is maintained.

The unit-sphere native and Python Cases are also selectable. Both run the
complete Mesh-to-CFD smoke chain, including the 3-D OneFLOW runtime. Because
the closed sphere has only one `wall` Physical Group (no inlet, outlet or
farfield), this is a runtime-pipeline smoke test, not an external-flow model or
a CFD accuracy benchmark. The staged flat-plate template supplies common
solver parameters only; the mesh dimension is overridden to 3.

From the repository root, run one selected Case after building OneFLOW and
making Gmsh available on PATH. For the native unit sphere:

    python -m tools.oneflow_gui.mesh_validation --oneflow-executable build/bin/Debug/OneFLOW.exe --gmsh gmsh --mesh-case-id stl_to_gmsh_sphere_automatic_tetrahedra --max-steps 2

For the Python reference, use `stl_to_gmsh_python_sphere_automatic_tetrahedra`
as the `--mesh-case-id`. For the prolate spheroid, use either
`stl_to_gmsh_prolate_spheroid_automatic_tetrahedra` or
`stl_to_gmsh_python_prolate_spheroid_automatic_tetrahedra`. These commands
execute the selected mesh through the complete Mesh-to-CFD smoke chain.

The default validation workspace is under
`%USERPROFILE%/.oneflow/mesh-validation/` on Windows or
`$HOME/.oneflow/mesh-validation/` on Linux. Each stage has its own directory,
logs and outputs. The source Cases are copied to the external workspace and
are never used as runtime directories.

The dedicated `grid/stl_to_gmsh/pipeline/` tree and its fixed Square workflow
have been retired. The STL Case runner below now runs only the
standalone STL-to-Gmsh producer Cases.

## Development tests

    python -m unittest grid.tests.python.stl_to_gmsh.test_stl_to_gmsh_cases
    python tools/mesh_cases.py validate stl_to_gmsh_square
    python tools/mesh_cases.py validate stl_to_gmsh_square_two_triangles
    python tools/mesh_cases.py validate stl_to_gmsh_square_2x2
    python tools/mesh_cases.py validate stl_to_gmsh_two_square_cells
    python tools/mesh_cases.py validate 2d/circle_external_flow/cpp

When Gmsh is available, the Python test suite executes the circle reference
Case, all three cube reference Cases, the prolate-spheroid reference Case and
the internal-flow cylinder reference Case in
temporary directories, validating their MSH 4.1 element families, fixed or
minimum topology contracts, coordinate bounds and Physical Groups. The native C++ integration suite separately exercises the
single-hex, six-tetrahedron, and automatic-tetrahedron cube strategies. The native C++ circle integration test is enabled
by `ONEFLOW_ENABLE_GMSH_INTEGRATION_TESTS` and is exercised by the
`OneFLOW External Tool Integration` workflow.

Tests use temporary output directories. The structured native Cases support
one connected, axis-aligned rectangular XY boundary and QUAD_4 layouts. The
native circle external-flow Case additionally supports one closed planar STL
body loop inside a rectangular outer domain, producing a triangular fluid
region with a hole. The three 3d/cube/* pairs remain the structured 3-D volume prototypes and are
restricted to the same axis-aligned cube STL. They progress from one HEXA_8
cell, to a fixed six-TETRA_4 decomposition, to an automatic tetrahedral mesh.
The separate 6:1 prolate-spheroid pair exercises the generic closed, connected,
consistently oriented genus-zero STL volume-meshing path. This does not claim
support for disconnected shells, non-zero genus surfaces, self-intersections,
or arbitrary invalid STL input.

## Retired legacy entry point

The former root-level `script/run.py` workflow was removed because it wrote generated files into the source Case tree. Use the registered `2d/circle_external_flow/cpp` and `2d/circle_external_flow/python` Cases through the isolated Case runner instead.

## 3D cube volume-mesh progression

The `3d/cube/` family holds one axis-aligned cube STL and six semantic boundary
patches constant while changing only the volume-mesh strategy. Start with
`single_hex` for one HEXA_8 cell, continue to `six_tetrahedra` for an exact
six-cell tetrahedral decomposition, and then use `automatic_tetrahedra` to
let Gmsh choose an unstructured tetrahedral volume mesh. See
[`3d/cube/README.md`](cube/README.md) for each Case's acceptance contract.


## 3D unit-sphere volume-mesh progression

The `3d/sphere/automatic_tetrahedra` native and Python Cases use one byte-identical,
watertight unit-sphere STL (center `(0, 0, 0)`, radius `1`). Gmsh classifies
the curved STL facets, reconstructs parametrizable patches, creates one closed
surface loop and meshes its interior with first-order tetrahedra. The runner
reuses this exact STL at coarse (0.40), medium (0.20) and fine (0.10) target
characteristic lengths for both implementations. These are mesh-size controls,
not requested element counts. The sphere surface is one `wall` Physical
Group and the interior is `fluid_domain`. See the 3D volume cases in
`../../docs/source/stl_to_gmsh_3d.rst` for the geometric contract and validation scope.


## 3D internal-flow cylinder asset

The `3d/cylinder_internal_flow/automatic_tetrahedra` Case pair is the
canonical internal-flow cylinder entry in the grid database. Native C++ and
Python reference Cases use the same analytic definition, radius/length and
48 x 16 surface discretization. Each representation carries its own checked-in
STL serialization and companion `.geometry.json` manifest.
The analytic definition is radius 1.0 m, length 4.0 m, axis X, centered at the
origin. Its region-aware Gmsh adapter creates separate `wall`, `inflow` and
`outflow` Physical Surfaces and the `fluid_domain` Physical Volume. The end
caps close the geometric shell; they are not wall patches.

The Mesh Generation preset references this same checked-in STL. The Geometry
Library may still generate custom cylinder variants into the user's
`.oneflow` workspace, but the canonical asset is never generated into or
modified in the source checkout. Both implementations are registered in the normal Case catalog and staged
externally by `run_cases.py`. The C++ implementation validates the closed
surface topology and manifest patch partition before constructing the OCC
cylinder; its validation is intentionally independent of the Python adapter.


## 3D curved pipe-bend geometry

The first curved-pipe asset is registered as `geometry.pipe_bend.90_degree`
in `grid/geometry/catalog.json`. Its independent Python and C++17 geometry
generators live under `3d/pipe_bend/90_degree/geometry/`. The Python generator
writes a checksum-bound closed STL, wall/inflow/outflow patch STL files and a
`region_manifest.json`; the existing generic region-aware volume-mesh route
consumes those files without pipe-bend-specific geometry logic. The GUI can
import the generated STL into a user mesh project and open the Gmsh workflow.
The C++17 generator writes the same named wall/inflow/outflow patch
STLs and checksum-bound `region_manifest.json` used by the Python path. The
native Case adapter compiles that generator in the external workspace and
hands its assets to the shared region-aware Gmsh volume mesher. See
`../../docs/source/stl_pipe_bend.rst` and the geometry README for equations,
defaults, commands and verification scope.
