6:1 Prolate Spheroid STL to Gmsh
=================================

Purpose
-------

The 6:1 prolate spheroid is a representative elongated, closed 3-D geometry
for exercising generic STL surface reconstruction and tetrahedral volume
meshing. This is a mesh-generation example and pipeline smoke case, not a
complete aerodynamic benchmark.

Geometry
--------

The spheroid is centered at the origin with semi-axes 3, 0.5 and 0.5:

.. math::

   x^2/9 + y^2/0.25 + z^2/0.25 = 1

Its full length is 6, maximum diameter is 1, aspect ratio is 6:1, and analytic
volume is pi. The checked-in STL scales the triangulated unit-sphere surface
by (3, 0.5, 0.5). Topology and facet count remain constant while proportions
change. The surface construction, affine transform, units and current
generator-provenance gap are recorded in :doc:`stl_geometry_generation`.

Case registration
-----------------

The native C++ Case is registered as
stl_to_gmsh_prolate_spheroid_automatic_tetrahedra and the Python reference Case
as stl_to_gmsh_python_prolate_spheroid_automatic_tetrahedra. Both use the same
byte-identical STL and equivalent layout contracts. The native Case appears
under Run Development Cases; the independent Python implementation appears
under Python reference cases. Both are available in the 3-D Mesh-to-CFD
Validation selection tree.

Mesh-to-CFD Validation
----------------------

The Case uses the existing STL -> Gmsh -> CGNS -> OneFLOW -> CFD smoke graph.
The default short 3-D runtime check verifies mesh handoff and solver execution.
The STL has only a wall Physical Group; this does not define an external-flow
problem or a physical accuracy, grid-independence or drag-coefficient test.
