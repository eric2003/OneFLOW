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
