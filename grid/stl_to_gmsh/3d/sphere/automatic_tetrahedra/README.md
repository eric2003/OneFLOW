# STL to Gmsh: Unit Sphere with Automatic TETRA_4 Meshing

This Case uses one closed, outward-oriented unit-sphere STL surface. Gmsh
classifies the STL facets, reconstructs parametrizable surface patches, forms
one closed surface loop and generates a tetrahedral volume mesh.

## Generate a volume mesh from an STL

Use the direct command; no OneFLOW Case directory is required:

```powershell
python grid/stl_to_gmsh/generate_volume_mesh.py --stl path/to/sphere.stl --output $HOME/.oneflow/mesh-generation/sphere.msh --mesh-size 0.20
```

Set a size interval instead of one uniform target:

```powershell
python grid/stl_to_gmsh/generate_volume_mesh.py --stl path/to/sphere.stl --output $HOME/.oneflow/mesh-generation/sphere.msh --mesh-size-min 0.10 --mesh-size-max 0.20
```

For coarse / medium / fine meshes, reuse the same STL and change only
`--mesh-size` to `0.40`, `0.20`, and `0.10`. Each run writes the
`.msh` mesh, a matching `.geo` file, and a `.report.json` validation
report. See [the concise command reference](../../../../../doc/source/stl_to_gmsh_commands.rst).

The sphere is centered at the origin with radius 1. Its input geometry is
fixed across all mesh-density variants. Only `mesh_size` changes:
coarse 0.40, medium 0.20 and fine 0.10. These values are target characteristic
lengths, not exact element counts.

The native C++ and independent Python reference Cases use byte-identical
STL inputs and the same layout contract. Generated files are written only to
the external staging directory by `grid/stl_to_gmsh/run_cases.py`.

The acceptance contract requires one connected, closed, consistently oriented
genus-zero STL surface, non-empty first-order TETRA_4 volume cells, matching
TRI_3 exterior facets, no QUAD_4/HEXA_8 cells, and the `wall` and
`fluid_domain` Physical Groups. Exact node and element counts are deliberately
not fixed because they depend on Gmsh version and meshing algorithm.

For the complete generation method, commands, STL input contract and validation limits, see [`stl_to_gmsh_sphere`](../../../../../doc/source/stl_to_gmsh_sphere.rst) and the Mesh Knowledge entry `mesh/stl-unit-sphere-to-gmsh`.
