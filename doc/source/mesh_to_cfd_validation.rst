Mesh-to-CFD Validation
======================

Mesh-to-CFD Validation is a dedicated mode in the OneFLOW Tests dialog. It is
separate from **Run Development Cases**, which continues to execute individual
preconfigured Cases.

The initial conversion route is:

.. code-block:: text

   selected STL-to-Gmsh Case -> Gmsh MSH -> Gmsh-to-CGNS Case
       -> CGNS-to-OneFLOW Case -> selected CFD template

The first CFD consumer is the existing ``test/cfd/plateuns2dslau2`` Laminar
Flat Plate Case. Each source Case is copied into a unique directory under
``~/.oneflow/mesh-validation`` (or the equivalent user-home directory on
Windows). Runtime parameter changes are applied only to those copies. The
runner derives the mesh artifact path from the selected Case's ``grid.txt``;
users do not edit workflow JSON or manually copy intermediate files.

The default solver smoke test advances two steps. A square mesh with the
flat-plate template tests conversion, native-grid ingestion, initialization,
and short time advancement; it is not a physical flat-plate validation.
Each stage writes a log, and the run writes ``mesh-validation-report.json``.
Additional CFD templates, including a turbulent flat plate, can be registered
as this mode evolves.
