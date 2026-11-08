# STL Shared Vertices within Weld Tolerance - native OneFLOW Case

This native C++ Grid Case repeats one shared facet vertex with a small
coordinate perturbation. The configured `stlWeldTolerance` is 1.0e-10,
while the repeated X coordinate differs by 2.0e-11.

Run from the repository root:

    OneFLOW.exe r grid/stl/weld_tolerance

The conversion should quantize the near-coincident coordinates to the same
vertex key and write two TRI_3 cells with four nodes.

PASS establishes the current configured tolerance-based vertex identity
behavior on this controlled fixture. It does not establish an optimal
tolerance for arbitrary scales, robust spatial clustering for every point
distribution, geometric repair, manifoldness or watertightness.
