# STL to Gmsh: Unit Sphere with Automatic TETRA_4 Meshing

This Case uses one closed, outward-oriented unit-sphere STL surface. Gmsh
classifies the STL facets, reconstructs parametrizable surface patches, forms
one closed surface loop and generates a tetrahedral volume mesh.

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

For the complete generation method, commands, STL input contract and validation limits, see [`stl_to_gmsh_sphere`](../../../../../../doc/source/stl_to_gmsh_sphere.rst) and the Mesh Knowledge entry `mesh/stl-unit-sphere-to-gmsh`.
