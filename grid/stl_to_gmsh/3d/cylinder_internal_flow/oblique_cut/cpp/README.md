# Cylinder Internal Flow — Native C++ Oblique Cut

This is a native OneFLOW grid-conversion case. It contains only case inputs and
configuration; the STL cutting algorithm is implemented in
`codes/geometry/stl`.

The case enables the generic cutter through `stlSurfaceCutEnabled = 1` in
`script/grid.txt`. Its shared cutting-plane configuration is
`script/cutting.json` (included in this self-contained native case workspace). The source STL is a blind triangulated shell; the
cutter does not read analytic cylinder metadata.

Run this case with the normal OneFLOW executable and case runner. The grid
conversion reads the STL, clips it against the configured inlet and outlet
planes, caps both loops, validates the closed shell, preserves the wall/inflow/
outflow surface regions, and invokes Gmsh for tetrahedral volume meshing.

Generated files are written under this case's `results/` directory. No
standalone C++ executable or algorithm source is stored in this case folder.
