# STL Vertex Welding across a Quantization Boundary - native OneFLOW Case

This native C++ Grid Case repeats a shared facet vertex with a coordinate
perturbation of 5.1e-11 while the configured stlWeldTolerance is 1.0e-10.

The two X coordinates lie on opposite sides of a rounding-cell boundary:
coordinate-wise rounding would assign different keys even though their
Euclidean distance is strictly smaller than the configured tolerance.

Run from the repository root:

    OneFLOW.exe r grid/stl_import/weld_tolerance

The conversion must use the exact spatial tolerance predicate and write two
TRI_3 cells with four unique nodes.

PASS exercises the production STL-to-CGNS conversion path for a tolerance
boundary case. It does not establish an optimal tolerance for arbitrary
scales, transitive clustering, geometric repair, manifoldness or watertightness.
