# STL ASCII Basic - native OneFLOW Case

This is a normal OneFLOW grid Case, not a Python utility case or a GoogleTest.
Run it through the standard application entry point:

    OneFLOW.exe r grid/stl/ascii_basic

The Case configuration in script/grid.txt selects the native C++ STL reader
through sourceGridType = "stl" and writes an unstructured surface CGNS file
to results/model.cgns.

stlWeldTolerance controls coordinate-based vertex welding. The bundled
two-triangle square should produce four mesh nodes and two TRI_3 cells.
The conversion writes results/model.stl-read-report.json.

This first Case converts a triangle surface. It does not generate a volume
mesh, infer units, or invent physical boundary groups.
