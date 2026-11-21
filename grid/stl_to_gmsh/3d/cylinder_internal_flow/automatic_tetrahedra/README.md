# Internal-flow cylinder: automatic tetrahedral volume mesh

This is the canonical 3D internal-flow cylinder geometry asset in the OneFLOW
grid database. It follows the same checked-in input convention as the cube and
unit-sphere STL-to-Gmsh Cases.

## Analytic definition

- Axis: X
- Center: origin
- Radius: 1.0 m
- Length: 4.0 m
- Axial extent: -2.0 m to +2.0 m
- Circumferential divisions: 48
- Axial divisions: 16
- Surface patches: `wall`, `inflow`, `outflow`

The STL is a closed geometric shell. The two end caps close the shell for
volume reconstruction; they are not wall boundaries. The region-aware Gmsh
OCC construction creates independent `wall`, `inflow`, `outflow`
Physical Surfaces and one `fluid_domain` Physical Volume.

## Inputs and outputs

Both native C++ and Python reference Cases are provided:

- `cpp/input/model.stl`: ASCII STL representation for the native C++ Case.
- `cpp/input/model.geometry.json`: analytic parameters, patch indices and
  SHA-256 binding for the ASCII representation.
- `cpp/script/layout.json`: native C++ volume-mesh strategy and density.
- `python/input/model.stl`: binary STL representation for the Python adapter.
- `python/input/model.geometry.json`: corresponding checksum-bound manifest.
- `python/script/layout.json` and `python/expected.json`: reference strategy
  and topology/Physical Group acceptance contract.

The two representations have the same analytic dimensions, triangle ordering
and patch partition; their byte checksums are representation-specific. Both
implementations are registered as independent Cases and can be run through
the repository runner. They stage under the user's
`.oneflow/grid-development/stl-to-gmsh` workspace; generated `.geo`, `.msh`,
reports and logs never enter the source tree.

The Geometry Library / Generate STL Geometry editor uses the same analytic
asset ID and default dimensions. Its user-selected parameterized outputs are
separate from this immutable standard database fixture.


## Stage A: blind STL surface inspection

Before using the analytic manifest or any patch indices, inspect the Python
Case's `input/model.stl` as an unknown triangle surface:

```powershell
python grid/stl_import/inspect.py --input grid/stl_to_gmsh/3d/cylinder_internal_flow/automatic_tetrahedra/python/input/model.stl --output "$env:USERPROFILE/.oneflow/grid-development/stl-inspection/cylinder_internal_flow/stage_a_report.json"
```

This Stage A command reads only `model.stl`. It must not read
`model.geometry.json`, use the `wall/inflow/outflow` facet lists, invoke
Gmsh, or generate a volume mesh. For the current canonical binary STL, the
expected topology is 1,632 triangles, one closed connected component, zero
boundary edges, zero non-manifold edges, Euler characteristic 2 and genus 0.
The report also records the signed surface volume and orientation. The
regression test enforces this blind inspection contract.

Stage B is intentionally separate: users will specify two cutting planes,
then the cut surface will be passed through Stage A again. The expected result
is two simple open boundary loops; their inlet/outlet meanings come from the
user's plane definitions, not from the original STL.


## Stage B: analytic cutting before generic STL cutting

The first cut is made from an explicit analytic cylinder definition:
axis X, centerline y=z=0, radius 1, extent -2 to +2. It is not reconstructed
from model.geometry.json and does not consume facet patch indices.

For planes x=-1.5 and x=+1.5, the equation y^2+z^2=R^2 gives two circular
intersection curves. Retain the lateral surface parameterized by
(x, R cos(theta), R sin(theta)) between those planes. Omit the original end
caps and do not create planar caps. The result is one open surface component
with exactly two boundary loops.

Run grid/stl_import/analytic_cylinder.py and then re-run the generic
grid/stl_import/inspect.py against the generated STL. Both outputs must be
under .oneflow/grid-development/stl-cutting/cylinder_internal_flow. The
default tessellation is 48 circumferential by 12 axial divisions.

This stage establishes the analytic reference only. It does not yet clip
arbitrary input STL triangles. Generic plane/triangle clipping on an unknown
curved pipe is the next stage, after the analytic boundary-loop checks pass.


## GUI cut-to-CFD workflow

Open **OneFLOW Workbench → Mesh Generation & Conversion → Cylinder Cutting & CFD**. The editor starts from the Geometry Library cylinder asset and exposes radius, axial extent, two cut-plane positions, surface divisions and Gmsh target size. The command preview shows the exact Python/Gmsh invocation.

Generation produces the open cut-surface STL, its Stage A report, the Gmsh geometry script and the MSH4 fluid-volume mesh, plus a combined generation report. Stage A must verify two open loops before Gmsh volume meshing is attempted.

The Gmsh volume is an analytic truncated cylinder, not a volume reconstructed from the open cut STL. It has explicit `wall`, `inflow`, `outflow` and `fluid_domain` Physical Groups. After mesh generation, **Run OneFLOW...** opens Mesh-to-CFD Validation with this MSH selected. That workflow performs MSH staging, CGNS conversion, OneFLOW grid ingestion and the registered short CFD smoke test. The smoke test does not claim cylinder-flow numerical accuracy; a dedicated internal-flow CFD Case remains separate work.
