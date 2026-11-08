OneFLOW Native STL Surface I/O
===============================

Purpose and status
------------------

This page defines the target architecture and verification contract for native
STL surface-file support in the C++ codes/geometry module.

**Status: proposed.** The current Python implementation in
grid/stl/stl_to_gmsh.py is a separate planar-surface conversion tool. Its tests
do not exercise OneFLOW's C++ geometry code. The native reader and writer
described here must be implemented and tested before the Mesh Contract
Workbench can report native STL I/O as implemented or verified.

STL's role in the geometry architecture
----------------------------------------

STL is a triangle-surface interchange format. It is not a volume mesh format,
and it does not carry the semantic information required to directly configure
a CFD solver. The ownership boundary is:

.. code-block:: text

   ASCII STL / Binary STL
             |
             v
      C++ StlMeshReader
             |
             v
        StlMesh value
             |
             +----> C++ StlMeshWriter ----> ASCII / Binary STL
             |
             v
      explicit surface validation
             |
             v
      explicit topology construction
             |
             v
      Geometry / Mesh Core adapters
             |
             +----> Gmsh / CGNS / OneFLOW mesh routes

The format adapter owns byte/text decoding, STL record validation and
serialization. It does not own vertex welding, manifold repair, boundary
condition assignment, volume meshing or solver mesh initialization. Those
operations belong to later, explicit services with their own contracts.

Proposed C++ module boundary
----------------------------

The implementation belongs beside the existing geometry format adapters:

.. code-block:: text

   codes/geometry/
       include/
           StlMeshIO.h
       src/
           StlMeshIO.cpp

The public API should be value-based and independent of Gmsh, CGNS, Qt and
Python. A first version can expose:

* StlFormat with ASCII, Binary and Auto values;
* StlTriangle: stored normal, three ordered vertices and a 16-bit binary
  attribute value;
* StlMesh: source header and ordered facet records;
* StlMeshReader::Read(path, format);
* StlMeshWriter::WriteASCII(path, mesh) and WriteBinary(path, mesh).

Exact names may be adjusted to match existing OneFLOW naming conventions, but
the ownership and data-preservation rules below are architectural invariants.
The reader should return a value or throw a documented parse/I/O exception;
it must not leave partially initialized global mesh state.

STL data semantics
------------------

ASCII STL
~~~~~~~~~

The conventional facet grammar is:

.. code-block:: text

   solid name
     facet normal nx ny nz
       outer loop
         vertex x1 y1 z1
         vertex x2 y2 z2
         vertex x3 y3 z3
       endloop
     endfacet
   endsolid name

The parser should accept ordinary whitespace and line-ending variation,
including LF and CRLF. Keywords should follow a documented case policy.
Numeric values must parse completely and be finite. The parser should reject
incomplete facets, missing terminators, malformed numbers and unsupported
facet structures with diagnostics that identify the file and parsing context.

ASCII STL has no standard per-facet attribute field. Writing a mesh that
contains nonzero binary attributes to ASCII is therefore lossy; the writer or
conversion report must state this rather than silently implying a lossless
round trip.

Binary STL
~~~~~~~~~~

The conventional binary layout is:

.. list-table::
   :header-rows: 1
   :widths: 20 20 60

   * - Offset
     - Size
     - Meaning
   * - 0
     - 80 bytes
     - Header (opaque bytes; not necessarily a null-terminated string)
   * - 80
     - 4 bytes
     - Unsigned triangle count, little-endian
   * - 84
     - 50 bytes per facet
     - 12 float32 values followed by one uint16 attribute value

Each facet contains three float32 normal components, nine float32 vertex
components and one uint16 attribute. The expected conventional file size is
84 + 50 * triangle_count. Implementations must check integer overflow before
computing the expected size, validate truncation, and use explicit little-endian
decoding/encoding. Do not rely on compiler struct packing or host endianness.

The 80-byte header is opaque. It may begin with "solid"; therefore, checking
the first five bytes is not sufficient to distinguish ASCII from binary.
Automatic detection should combine file length and structural parsing, and
must report ambiguous or malformed inputs instead of guessing silently.
Any policy for trailing bytes should be explicit and tested.

Facet identity, normals and topology
------------------------------------

STL repeats three vertex coordinates in every facet. It does not provide
authoritative shared node IDs. The reader must preserve facet order and vertex
order and must not weld coordinates as a side effect of parsing.

The stored normal and vertex winding are source records. A normal may be zero,
non-unit or inconsistent with the winding. Reading should preserve the source
normal; validation may report its magnitude and its agreement with the
geometric normal. Automatic normal recomputation, winding reversal and
coordinate welding are repair/transform operations and must never happen
silently inside the reader.

For vertices x1, x2, x3, the geometric area is

.. math::

   A = \frac{1}{2}
       \left\|
       (\mathbf{x}_2-\mathbf{x}_1) \times
       (\mathbf{x}_3-\mathbf{x}_1)
       \right\|.

A zero or tolerance-degenerate area is a geometric diagnostic. Whether such a
facet is rejected at parse time or retained with a validation finding must be
an explicit API policy; the recommended separation is to parse finite records
first and let a surface validator classify degeneracy.

STL does not standardize units, coordinate reference systems, material
semantics, named boundary patches, physical groups, volume elements or solver
boundary conditions. Those must come from external project metadata or later
mapping contracts. A file being syntactically valid does not imply a closed,
manifold or CFD-ready surface.

Reader/writer acceptance criteria
---------------------------------

.. list-table::
   :header-rows: 1
   :widths: 24 48 28

   * - Area
     - Required behavior
     - Evidence
   * - ASCII reader
     - Basic facet, whitespace, CRLF, numeric notation, malformed and truncated records
     - C++ GoogleTests
   * - Binary reader
     - Header, count, little-endian values, attributes, solid-prefixed header, truncation and count/length mismatch
     - C++ GoogleTests
   * - Format detection
     - ASCII/binary discrimination based on structural evidence, including ambiguous headers
     - C++ GoogleTests
   * - ASCII writer
     - Stable precision, valid grammar, finite-value policy
     - Parse-after-write test
   * - Binary writer
     - Exact header/count/facet layout, little-endian output, attribute preservation
     - Byte-level assertions and parse-after-write test
   * - Round trip
     - Compare facet count, order, coordinates, normals and all representable attributes
     - ASCII and binary round-trip tests
   * - Error handling
     - File-open failures, invalid numbers, truncated records, unsupported structure, overflow and output failures
     - Exception/diagnostic tests
   * - Cross-platform
     - No host-endian or Windows-only path assumptions
     - Windows and Linux CI

The initial tests should be native C++ tests linked to OneFLOW_GeometryCore.
A fixture may be shared with Python, but a Python test invoking
stl_to_gmsh.py is not evidence for this C++ contract. Integration tests that
build shared topology or convert STL to a OneFLOW mesh are a later milestone
and must be reported separately from parser tests.

Relationship to existing Python tools
-------------------------------------

The existing grid/stl Python utilities remain useful for planar boundary
extraction and Gmsh geometry generation. They have a different responsibility:

* Python utility: planar boundary extraction, sidecar-driven domain setup and
  Gmsh generation.
* C++ Geometry adapter: native STL file I/O and source-facet preservation.
* Later topology service: explicit vertex welding, connected components,
  edge incidence, manifold checks and surface orientation diagnostics.
* Later conversion adapter: validated surface topology into a mesher or
  OneFLOW Mesh Core representation.

The eight Python STL development cases currently registered in the development
case catalog exercise the Python utility only. They must not be described as
native C++ STL reader tests. Native C++ tests should be added to CTest and
exposed through the appropriate development-test inventory after the C++
implementation exists.

Mesh Contract Workbench entry
-----------------------------

The machine-readable contract is
tools/verification/contracts/mesh/stl-native-io.json, with ID
mesh/stl-native-io. It is grouped under Mesh / Interchange Formats and
currently has status proposed. The Overview, Principle, Reference and
Architecture tabs describe the target. Implementation source paths and
verification evidence remain empty until corresponding C++ files and tests
are committed. This prevents a design declaration from being mistaken for
working code.

References
----------

* Library of Congress, STL (STereoLithography) File Format:
  https://www.loc.gov/preservation/digital/formats/fdd/fdd000504.shtml
* Gmsh Reference Manual, file formats and geometry import:
  https://gmsh.info/doc/texinfo/gmsh.html
* OneFLOW mesh_architecture, gmsh_interface and stl_to_gmsh_2d pages.
