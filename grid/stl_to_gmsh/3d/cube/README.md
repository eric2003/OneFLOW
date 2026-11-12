# 3D cube STL-to-Gmsh progression

All three Cases use the same closed, axis-aligned cube STL. They differ only
in the volume-mesh strategy, so each stage isolates one topology question.

| Case | Gmsh strategy | Expected volume topology | Purpose |
|---|---|---|---|
| [single_hex](single_hex/README.md) | Transfinite + recombination | 1 HEXA_8 | Minimum hexahedral volume |
| [six_tetrahedra](six_tetrahedra/README.md) | Transfinite, no recombination | 6 TETRA_4 | Fixed tetrahedral subdivision |
| [automatic_tetrahedra](automatic_tetrahedra/README.md) | Automatic unstructured volume mesh | Gmsh-generated TETRA_4 count | Automatic tetrahedral meshing |

Each Case has native C++ and independent Python reference implementations.
Both consume byte-identical STL geometry and equivalent layout contracts.
They preserve the six boundary Physical Groups
(`xmin_inflow`, `xmax_outflow`, `ymin_wall`, `ymax_farfield`,
`zmin_symmetry`, `zmax_symmetry`) and the `fluid_domain` volume group.

The six-tetrahedron Case has a fixed acceptance contract: eight corner nodes,
twelve boundary TRI_3 facets and exactly six TETRA_4 volume cells. The
automatic Case does not prescribe node or cell counts; it requires a
non-empty tetrahedral volume mesh and triangular boundary facets, and rejects
quadrilateral or hexahedral output.

Run the registered native and Python Cases from the repository root:

```sh
python grid/stl_to_gmsh/develop_test/run_cases.py --oneflow-executable build/bin/Debug/OneFLOW.exe --gmsh gmsh
```

The runner stages all Cases outside the source tree under the user's
`.oneflow/mesh-layout-test/grid/stl_to_gmsh` directory.
