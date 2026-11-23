STL surface diagnostics
=======================

Purpose and boundaries
----------------------

Stage A provides general STL import and non-destructive surface diagnostics.
It is independent of Gmsh and does not prepare inlet/outlet boundaries. The
parser preserves triangle records; diagnostics construct temporary
shared-vertex topology using a documented tolerance. The source STL is never
modified.

.. code-block:: console

   python grid/stl_import/inspect.py --input model.stl --output stl_report.json

The report includes representation, triangle count, coordinate bounds,
bounding-box diagonal, surface area, welded vertex and edge counts,
degenerate and duplicate facets, boundary and non-manifold edges, orientation
conflicts, connected components, and closedness. For closed, consistently
oriented components it also reports Euler characteristic, genus when
applicable, signed volume, and outward winding.

Tolerance and units
-------------------

STL does not store a reliable unit declaration. Coordinates are reported as
read and units are marked unknown. The default weld tolerance is 1e-9 times
the bounding-box diagonal with a 1e-12 numerical floor. An explicit absolute
tolerance in model units can be supplied with --weld-tolerance. The selected
tolerance is recorded because it affects vertex identity and topology.

Interpretation
--------------

Readable means the file representation could be parsed; it does not mean the
surface is suitable for volume meshing.

* valid: no detected degenerate or duplicate facets, non-manifold edges, or
  inconsistent edge orientation; one connected component.
* warning: readable but open or multi-component. Open boundaries can be
  intentional, including an uncapped tube surface.
* invalid: degenerate/duplicate facets, non-manifold edges, or inconsistent
  orientation were detected.

volume_mesh_ready is stricter: one closed, edge-manifold, consistently
oriented component with positive signed volume. Genus zero is not required;
genus is reported rather than assumed.

Limitations
-----------

This stage does not test triangle-triangle self-intersections, repair
intersections, fill holes, reorient facets, infer units, classify cuts as
inlets/outlets, assign CFD boundary conditions, or generate a volume mesh.
A clean report is not proof that the surface is self-intersection-free.

For an internal-flow cylinder, the next stage can cut the surface at
user-specified planes. Stage A should report the resulting open boundary
loops; it must not silently close them or label them inlet/outlet.

Implementation
--------------

Implementation: grid/stl_import/reader.py and grid/stl_import/validation.py.
Python reference tests: grid/tests/python/stl_import/. Native C++
StlMeshReader continues to own native ASCII/Binary parsing; this report does
not claim native C++ topology-validation parity.
