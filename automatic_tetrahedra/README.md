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

- `python/input/model.stl`: canonical, version-controlled surface geometry.
- `python/input/model.geometry.json`: asset identity, analytic parameters,
  patch triangle indices and SHA-256 binding to the STL.
- `python/script/layout.json`: mesh-generation strategy and density.
- `python/expected.json`: topology and Physical Group acceptance contract.

Run the Python reference Case through the repository runner. It stages the
Case under the user's `.oneflow/grid-development/stl-to-gmsh` workspace;
generated `.geo`, `.msh`, reports and logs never enter the source tree.

The Geometry Library / Generate STL Geometry editor uses the same analytic
asset ID and default dimensions. Its user-selected parameterized outputs are
separate from this immutable standard database fixture.
