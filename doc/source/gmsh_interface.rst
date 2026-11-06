Gmsh MSH Interface
==================

Purpose and implementation status
----------------------------------

OneFLOW reads the Gmsh native mesh file format (MSH) through a file-format
adapter. The adapter is intentionally independent of the Gmsh SDK and API:
the solver does not need to link against Gmsh, and users can generate meshes
with a separately installed Gmsh application.

The current implementation is an input-validation prototype connected to the
normal Case parameter path. It does not yet construct a canonical OneFLOW
mesh from the parsed coordinates and connectivity, and its current
.cgns-named output is a marked text placeholder, not a CGNS file.

Version policy
--------------

MSH 4.1 ASCII is the primary target. MSH 4.1 is the current 4.x file-format
revision documented by Gmsh. The legacy MSH 2.x ASCII reader is retained as a
compatibility path. The supported matrix is deliberately narrower than the
complete Gmsh format:

.. list-table::
   :header-rows: 1
   :widths: 18 20 62

   * - Input
     - Status
     - Notes
   * - MSH 4.1 ASCII
     - Prototype supported
     - Entity-block node and element sections; sparse node tags; parametric
       node coordinates; common element types 1-31.
   * - MSH 2.x ASCII
     - Prototype supported
     - Legacy line-oriented node and element records.
   * - MSH 4.0
     - Not supported
     - Its block layout differs from MSH 4.1 and requires a separate reader.
   * - Binary MSH
     - Not supported
     - Binary integer widths, endianness, and block payloads require a
       separate, explicitly tested implementation.
   * - MSH 4.1 element types outside 1-31
     - Not supported
     - Variable-size polygon/polyhedron and later high-order types need
       explicit connectivity rules.

MSH 4.1 structure and parser responsibilities
----------------------------------------------

A MSH file is organized as named sections. The parser first checks
$MeshFormat and rejects binary input before interpreting text records.
For MSH 4.1, the mesh payload is grouped by model entities:

* $Entities describes geometric points, curves, surfaces and volumes,
  including their physical-group memberships and boundary relationships.
* $Nodes starts with the number of entity blocks, total node count, and
  minimum/maximum node tags. Each block identifies its entity dimension and
  tag, whether nodes are parametric, and its node count. The block then stores
  all node tags followed by coordinates (and parametric coordinates where
  applicable).
* $Elements also starts with block and total counts. Each block identifies
  its entity dimension, entity tag, element type and count. Element records
  refer to nodes by node tag, not by the node's position in a dense array.

Node and element tags are positive identifiers and are not guaranteed to be
contiguous or ordered. A production reader must therefore build a tag-to-local
index map; it must never treat a Gmsh tag as a C++ vector offset. The current
prototype validates uniqueness and node references using tag sets. It parses
the coordinates and connectivity records but does not yet retain them in the
OneFLOW mesh model.

The parser also needs to distinguish three different concepts:

* Geometric entity: the Gmsh model entity (point, curve, surface or volume)
  on which mesh elements and nodes are classified.
* Physical group: a user-defined semantic group, often used to identify
  fluid regions, walls, inlets, outlets, interfaces or material regions.
* Element type: the interpolation/topology type and node ordering of an
  individual mesh element.

These identifiers are not interchangeable. In particular, a geometric
entity tag is not automatically a solver boundary-condition ID. A future
mapping layer must resolve physical-group names and dimensions explicitly,
then map them into OneFLOW region/boundary contracts.

Current validation contract
----------------------------

The MSH 4.1 ASCII prototype currently checks:

* valid file-format version and ASCII mode;
* non-empty node and element sections and consistent block totals;
* positive, unique node and element tags;
* supported element type and its declared node count;
* each element's node references resolve to a declared node tag;
* section terminators and truncated records.

The MSH 2.x path retains its existing basic record/count validation. It is a
compatibility reader, not yet a semantic mesh importer.

Not yet preserved
-----------------

The following information is intentionally not claimed as supported:

* $Entities bounding representation and entity adjacency;
* $PhysicalNames and entity-to-physical-group membership;
* periodic node relations;
* partitioned entities, ghost elements and partition ownership;
* high-order interpolation metadata beyond the currently recognized
  element connectivity counts;
* post-processing sections such as $NodeData, $ElementData and
  $ElementNodeData;
* units, coordinate reference conventions and application-specific metadata.

The placeholder summary explicitly records that geometric entities and
physical groups are not yet preserved. This prevents a successful parse from
being mistaken for a lossless conversion.

Target architecture
-------------------

The intended dependency direction is:

.. code-block:: text

   Gmsh MSH 4.1 / MSH 2.x
            |
            v
      Gmsh file reader
            |
            v
   Format-neutral mesh snapshot
     - coordinates and stable node IDs
     - element types and ordered connectivity
     - geometric entity classification
     - physical groups and names
     - periodic / partition metadata when supported
            |
            v
     OneFLOW mesh adapter
            |
            +----> native OneFLOW mesh
            |
            +----> CGNS writer / reader
            |
            +----> validation and inspection reports

The format-neutral snapshot is the critical boundary. Gmsh parsing should not
depend on solver globals, GUI widgets, or CGNS library types. Likewise, the
CGNS adapter should not know how MSH sections are tokenized. The same snapshot
can later support Gmsh-to-CGNS, CGNS-to-Gmsh, and direct Gmsh-to-OneFLOW
conversion without implementing every pair as an unrelated converter.

Conversion validation strategy
-------------------------------

File creation is not a sufficient acceptance criterion. The planned
verification levels are:

#. Syntax and structural validation — section counts, tag uniqueness,
   block totals, element arity and reference validity.
#. Geometry validation — coordinate finiteness, bounding box, dimension,
   duplicate nodes and scale-aware geometric tolerances.
#. Topology validation — cell connectivity, face incidence, orientation,
   boundary classification and connected components.
#. Semantic validation — physical-group names and dimensions mapped to
   named OneFLOW regions and boundary patches.
#. Round-trip validation — export to another format, reload it, then
   compare coordinates, element types, ordered connectivity and semantic
   region assignments. Comparisons must account for permitted node/element
   renumbering rather than requiring identical file ordering.

The GUI exposes this prototype only as a manual experimental Case check.
It is not part of CTest, Python unittest, the regular CFD/grid suites, or
GitHub Actions. The Case is copied to a temporary workspace before execution;
the repository's sample input remains unchanged.

References
----------

* `Gmsh Reference Manual, MSH file format <https://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format>`_
* `Gmsh Reference Manual, current release <https://gmsh.info/doc/texinfo/>`_
* Geuzaine, C. and Remacle, J.-F., *Gmsh: A three-dimensional finite element
  mesh generator with built-in pre- and post-processing facilities*,
  International Journal for Numerical Methods in Engineering, 79(11),
  1309-1331, 2009.
