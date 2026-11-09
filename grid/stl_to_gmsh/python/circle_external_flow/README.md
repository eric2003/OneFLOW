# 2D Circle External Flow (Python Reference)

This independent Python reference uses the same 32-segment circular STL,
outer-domain bounds, mesh-size controls and Physical Group names as the native
C++ Case at `grid/stl_to_gmsh/circle_external_flow`. Both construct the fluid
surface as the rectangular outer loop with the STL body loop as a hole.

Run from the repository root:

    python grid/stl_to_gmsh/python/run_case.py --case-directory grid/stl_to_gmsh/python/circle_external_flow

The unstructured mesh is checked by Physical Group names and minimum node,
boundary-line and triangle counts rather than exact counts.
