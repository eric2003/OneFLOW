# STL Non-planar Surface - native OneFLOW Case

This native C++ Grid Case contains two adjacent triangles that form a folded
surface. It deliberately has no single constant-Z plane.

Run from the repository root:

    OneFLOW.exe r grid/stl/non_planar_surface

The standard Grid task reads the STL, welds the shared edge vertices and
writes an unstructured triangle surface to CGNS. Expected report:
two triangles and four welded nodes.

PASS establishes that the native STL-to-CGNS surface conversion preserves
three-dimensional coordinates and does not impose the Python planar-boundary
extractor's single-plane restriction. It does not prove manifoldness,
watertightness, orientation consistency, surface repair, volume meshing,
units, physical groups or CFD boundary-condition assignment.
