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

Gmsh API integration and cross-platform policy
----------------------------------------------

The current MSH reader does not link to Gmsh. This remains the preferred
baseline for importing existing mesh files: it keeps OneFLOW builds
self-contained, avoids a runtime Gmsh dependency, and makes the file-format
contract directly testable.

If OneFLOW later needs to create geometry, generate meshes, query CAD entities
or invoke Gmsh algorithms in-process, the official Gmsh SDK can be integrated
as an optional dependency. The SDK is distributed for Windows, Linux and
macOS and contains the shared library, headers and language bindings. This
does not require separate platform-specific OneFLOW API implementations:
the same adapter source can call the same Gmsh API on each platform.

The platform-specific pieces are the SDK artifact and its build/runtime
properties:

* Windows: Gmsh DLL plus its import library and headers; the DLL must be
  discoverable at runtime.
* Linux: Gmsh shared object plus headers; the shared library and its
  dependent libraries must be discoverable by the dynamic loader.
* macOS: Gmsh dynamic library plus headers and the platform loader/rpath
  configuration.

Compiler ABI is a separate concern from operating-system portability. The
native Gmsh C++ API may require an ABI-compatible compiler toolchain. The
Gmsh SDK documents a gmsh.h_cwrap compatibility header that exposes the
C++-style API through the C API for toolchains that cannot use the native
C++ ABI. For a long-lived multi-compiler project such as OneFLOW, the
preferred integration boundary is therefore the stable C API (or the
official C-API-backed wrapper), hidden behind a small OneFLOW-owned adapter.
Do not expose Gmsh types in public OneFLOW geometry headers.

The future integration should follow these rules:

#. Keep Gmsh support optional and disabled by default. Reading an existing
   MSH file must continue to work without the Gmsh SDK.
#. Add a dedicated CMake option (for example
   ONEFLOW_ENABLE_GMSH_API) only when an actual API-backed feature is
   implemented. Resolve the SDK through an explicit GMSH_ROOT hint and
   imported target; do not hard-code Visual Studio, GCC or installation paths.
#. Compile one shared adapter implementation. Platform-specific CMake logic
   should only select the SDK library/import library and runtime search path.
#. Keep SDK version checks and API capability checks in the adapter, not in
   solver or GUI code.
#. Convert Gmsh-owned data immediately into the format-neutral OneFLOW mesh
   representation. No Gmsh handles, pointers, entity classes or containers
   should leak into solver-facing mesh types.
#. Test both configurations: API disabled (the default, no Gmsh installation
   required) and API enabled (with an SDK explicitly supplied). CI should
   never download or require the SDK unless that build variant is deliberately
   selected.

This leads to a clean separation:

.. code-block:: text

   Optional Gmsh SDK (Windows DLL / Linux SO / macOS dylib)
                         |
                 Gmsh API adapter
                 (one C++ source)
                         |
                Format-neutral mesh
                         |
             OneFLOW geometry / mesh core

The standalone gmsh executable is another valid integration option when
OneFLOW only needs to launch meshing scripts or exchange files. It avoids
linking Gmsh into the OneFLOW process, at the cost of process management and
file-based exchange. The ONELAB socket interface is an option for interactive
parameter exchange. These are separate deployment choices, not reasons to
duplicate OneFLOW's geometry implementation by platform.

References for the API and SDK
------------------------------

* `Gmsh API and SDK installation <https://gmsh.info/doc/texinfo/gmsh.html#Installing-and-running-Gmsh-on-your-computer>`_
* `Gmsh API reference <https://gmsh.info/doc/texinfo/gmsh.html#Gmsh-application-programming-interface>`_
* `Gmsh SDK ABI notes <https://github.com/gmsh-project/gmsh/blob/master/doc/SDK.txt>`_

References
----------

* `Gmsh Reference Manual, MSH file format <https://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format>`_
* `Gmsh Reference Manual, current release <https://gmsh.info/doc/texinfo/>`_
* Geuzaine, C. and Remacle, J.-F., *Gmsh: A three-dimensional finite element
  mesh generator with built-in pre- and post-processing facilities*,
  International Journal for Numerical Methods in Engineering, 79(11),
  1309-1331, 2009.
