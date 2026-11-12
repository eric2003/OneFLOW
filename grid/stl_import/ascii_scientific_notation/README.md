# STL ASCII Scientific Notation - native OneFLOW Case

This is a native C++ OneFLOW Grid Case. It exercises the production STL
reader and STL-to-CGNS conversion path using finite coordinates written in
scientific notation (for example, `1.0e+00` and `1.0e-10`).

Run from the repository root:

    OneFLOW.exe r grid/stl/ascii_scientific_notation

The Case intentionally has no `script/control.txt`. The runtime discovers
`dimension.txt`, `grid.txt`, and `simu.txt` in the script directory.
The two-triangle square should produce four welded nodes and two TRI_3 cells.
The generated CGNS mesh and STL read report are written below `results/`.

Passing proves that the native runtime route accepts this ASCII numeric
notation and builds the expected triangle surface. It does not prove volume
meshing, units, boundary-condition semantics, manifoldness, or CFD readiness.
