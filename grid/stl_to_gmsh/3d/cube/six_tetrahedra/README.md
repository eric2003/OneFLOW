# STL to Gmsh: Cube with Six TETRA_4 Cells

This Case holds the cube geometry fixed and changes only the volume topology
from one HEXA_8 cell to a deterministic six-tetrahedron subdivision.

The input STL is byte-identical to the `single_hex` Case. Gmsh receives the
same eight cube corners, six planar boundary surfaces and six boundary
Physical Groups. Curves, surfaces and volume are transfinite, but surface and
volume recombination are disabled. The acceptance contract is exactly eight
nodes, twelve boundary TRI_3 facets and six volume TETRA_4 cells.

This fixed-count Case is not the automatic tetrahedral mesher. It isolates
the conversion and connectivity path for a known, reproducible tetrahedral
decomposition. The generated MSH 4.1 mesh is checked for exact element
families and counts, and all six boundary Physical Groups must survive.
The tetrahedral mesh is also checked topologically: every exterior face must
have a matching TRI_3 element, every TRI_3 must belong to a tetrahedral face,
and no face may be shared by more than two tetrahedra. This check uses node
connectivity rather than a prescribed cube diagonal.

Native C++ and Python reference Cases share the same STL and layout. See
[the cube progression](../README.md) and
[Development Case Methodology](../../../../../doc/source/development_case_methodology.rst).
