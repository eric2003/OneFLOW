# STL Binary Header Starts with solid - native OneFLOW Case

This is a native C++ OneFLOW Grid Case. It exercises binary STL input through
the standard Grid task and STL-to-CGNS conversion path. The fixture contains
two triangles forming a square and should produce four welded nodes and two
TRI_3 cells.

Run from the repository root:

    OneFLOW.exe r grid/stl_import/binary_header_solid

The Case has no `script/control.txt`; runtime parameters are discovered from
`script/dimension.txt`, `script/grid.txt`, and `script/simu.txt`.
The `binary_header_solid` fixture is a conventional 184-byte binary STL containing an
80-byte header, a two-facet count, and two 50-byte facet records.
Its opaque 80-byte header begins with the ASCII characters `solid`; this specifically prevents format detection from relying on a five-byte prefix test.

Passing proves that the native C++ reader and runtime conversion path decode
this binary fixture and produce a surface CGNS output. It does not prove
volume meshing, unit interpretation, physical groups, boundary-condition
mapping, manifoldness, or watertightness.
