# Format-owned mesh fixtures

Store each external format's smallest representative input here, once.
Conversion routes must reference these fixtures instead of copying them.

During the initial migration, the canonical Gmsh MSH 4.1 square fixture
remains at `grid/conversion/gmsh_to_cgns/grid/square.msh`. It is not duplicated here.
A future move must update all CMake test paths and manual Case references in
the same functional change.

Each format note should state supported versions, ASCII/binary scope,
element families, ID conventions, metadata preserved and known limitations.
