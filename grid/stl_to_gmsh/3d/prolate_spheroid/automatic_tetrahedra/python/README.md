6:1 Prolate Spheroid - Python Reference
========================================

Purpose
-------

This Case exercises generic closed, connected, consistently oriented, genus-zero
STL-to-Gmsh volume meshing on an elongated body. It has full length 6, maximum
diameter 1 and aspect ratio 6:1, centered at the origin.

Geometry
--------

.. math::

   x^2/3^2 + y^2/0.5^2 + z^2/0.5^2 = 1

The semi-axes are (3, 0.5, 0.5), and the analytic volume is pi. The checked-in
STL scales the triangulated unit-sphere surface by (3, 0.5, 0.5), preserving
its connected genus-zero triangulation while changing the body proportions.
It is a piecewise-planar approximation, not an exact CAD surface.

Mesh contract
-------------

* One closed, connected, consistently oriented genus-zero surface.
* One tetrahedral volume with Physical Surface wall and Physical Volume fluid_domain.
* First-order TETRA_4 volume elements and TRI_3 exterior facets.
* Mesh size is a target characteristic length, not an exact element count.

Mesh-to-CFD scope
-----------------

This Case is registered for Mesh-to-CFD Validation and runs the same short
3-D runtime smoke chain as the sphere. Its only Physical Group is a wall; it
does not define inlet, outlet, freestream, angle of attack or a boundary-layer
mesh. This is a mesh handoff/runtime smoke test, not an external-flow accuracy
benchmark or a validated spheroid drag coefficient.
