OneFLOW Native STL Surface I/O
===============================

Purpose and status
------------------

This page defines the target architecture and verification contract for native
STL surface-file support in the C++ codes/geometry module.

**Status: first C++ I/O slice implemented; CI verification pending.** The current Python implementation in
grid/stl/stl_to_gmsh.py is a separate planar-surface conversion tool. Its tests
do not exercise OneFLOW's C++ geometry code. The native reader and writer
are now implemented in StlMeshIO.h and StlMeshIO.cpp. Native GoogleTests are
registered in tests/geometry/StlMeshIOTest.cpp; their execution status is
tracked separately by the Mesh Contract Workbench.

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

Development Case and verification matrix
----------------------------------------

The Development Tests catalog separates native OneFLOW Cases from Python
reference Cases. In the Tests dialog, the **Python examples** checkbox is an
execution-mode selector:

* Unchecked (default): run only native C++ OneFLOW Cases through the normal
  `OneFLOW.exe r <case-directory>` entry point.
* Checked: run only Python reference Cases through their declared Python
  script. No OneFLOW executable is required in this mode.
* The case selection is shared and retained when switching modes. Each Run
  executes only the selected cases belonging to the active mode.

This separation is essential: a Python parser passing says nothing about
whether `StlMeshReader`, the C++ Grid conversion dispatcher, or the CGNS
writer works. Conversely, a native runtime Case does not validate the
independent Python geometry helper.

Native C++ OneFLOW Development Cases
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 23 31 24 22

   * - Case
     - Purpose
     - What PASS establishes
     - Not covered
   * - `stl_ascii_basic`
     - Baseline ASCII two-triangle square through the production Grid task
     - C++ STL read, tolerance-based shared-vertex construction, TRI_3 surface and CGNS output are connected; expected four nodes and two triangles
     - Volume meshing, units, physical groups, BC assignment, manifold or watertight proof
   * - `stl_ascii_crlf`
     - Same geometry with CRLF line endings
     - Native C++ ASCII parsing is independent of Windows/Linux text line endings
     - Other encodings, arbitrary malformed grammar, geometry repair or volume generation
   * - `stl_ascii_scientific_notation`
     - Same geometry with exponent-form finite coordinates and normals
     - Native C++ numeric token parsing accepts scientific notation through the runtime conversion route
     - Extreme exponent/overflow boundaries, locale-specific decimal separators, geometric quality or solver suitability
   * - `stl_binary_basic`
     - Two-triangle square encoded as standard binary STL
     - Native C++ runtime decodes the 80-byte header, little-endian count and 50-byte facets and writes the surface CGNS output
     - All binary corruption cases, arbitrary vendor extensions, attribute semantic interpretation, volume meshing
   * - `stl_binary_header_solid`
     - Binary STL whose opaque header begins with the word `solid`
     - Native C++ runtime detects binary structure rather than relying on the header prefix
     - Every ambiguous/trailing-byte variant or all possible format-detection ambiguities
   * - `stl_non_planar_surface`
     - Two adjacent facets form a folded, non-planar 3D surface
     - Native C++ STL-to-CGNS conversion preserves Z coordinates without assuming a single planar XY boundary
     - Surface repair, self-intersection detection, manifoldness, volume meshing or solver boundary assignment
   * - `stl_weld_tolerance`
     - Adjacent facets repeat a shared vertex with a small coordinate perturbation
     - Configured `stlWeldTolerance` maps near-coincident coordinates to one node; expected output is four nodes and two triangles
     - Robust spatial clustering for arbitrary point distributions, tolerance sensitivity studies, and topology quality beyond vertex identity

Python reference Development Cases
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 23 31 24 22

   * - Case
     - Purpose
     - What PASS establishes
     - Not covered
   * - `stl_python_ascii_basic`
     - Baseline ASCII parser and planar-boundary extraction
     - Python utility reads the fixture and extracts the expected simple closed boundary
     - Native C++ parsing, CGNS writing, general 3D topology
   * - `stl_python_ascii_crlf`
     - ASCII line-ending compatibility
     - Python reference path handles CRLF input
     - Native C++ line-ending behavior or other encodings
   * - `stl_python_binary_basic`
     - Conventional 80-byte header and binary triangle record
     - Python utility recognizes and reads a structurally valid binary STL
     - C++ binary reader/writer, attribute round trip, large-file robustness
   * - `stl_python_binary_header_solid`
     - Binary header whose first bytes spell `solid`
     - Python format detection does not classify binary solely from the header prefix
     - Native C++ auto-detection (covered separately by C++ unit tests), all ambiguous files
   * - `stl_python_scientific_notation`
     - ASCII exponent notation
     - Python numeric parser accepts the fixture's finite scientific notation
     - Native C++ numeric parsing or arbitrary numeric edge cases
   * - `stl_python_degenerate_triangle`
     - Zero-area / collapsed triangle diagnostic in planar extraction
     - Python topology stage detects and reports this invalid geometry
     - Whether the C++ STL reader rejects it; parsing and geometric validation are separate contracts
   * - `stl_python_truncated_binary`
     - Incomplete binary STL payload
     - Python reader reports the expected truncated-input diagnostic rather than accepting partial data
     - C++ binary truncation behavior (covered by C++ unit test), all corrupt binary variants
   * - `stl_python_non_planar_surface`
     - Non-planar input to a planar-boundary workflow
     - Python utility refuses to treat unsupported non-planar geometry as a valid planar boundary
     - Surface remeshing, projection, repair, or general 3D surface validity

Native C++ GoogleTests
~~~~~~~~~~~~~~~~~~~~~~

The focused `StlMeshIOTest` suite checks the I/O API in isolation using
temporary fixtures. Its individual test purposes are:

* `ReadsAsciiFacetAndPreservesVertexOrder`: proves the ASCII reader returns
  the stored normal and ordered vertices; it does not prove CGNS conversion.
* `ReadsBinaryAndPreservesSolidHeaderAndAttribute`: proves binary layout
  decoding and preservation of the opaque header and 16-bit facet attribute;
  it does not prove topology construction or all vendor-specific attributes.
* `RejectsTruncatedBinaryFacet`: proves a size/count mismatch is rejected;
  it does not exhaust every malformed-file condition.
* `RejectsNonFiniteAsciiCoordinates`: proves non-finite numeric data is
  rejected at the I/O boundary; it does not assess finite but poor-quality
  geometry.
* `WritesBinaryAndRoundTripsAttributeAndFacetData`: proves binary output
  can be read back with facet data and attributes preserved; it does not
  guarantee byte-identical output for arbitrary source headers.
* `WritesAsciiAndRoundTripsFacetGeometry`: proves ASCII serialization
  preserves representable facet geometry; ASCII has no binary attribute field.
* `RejectsAsciiWriteWhenBinaryAttributesWouldBeLost`: proves the writer
  refuses a known lossy conversion; it does not define every downstream
  interchange policy.
* `RejectsBinaryWriterValuesOutsideFloat32Range`: proves binary output
  rejects values that cannot be represented as finite float32; it does not
  guarantee adequate precision for every engineering scale.

For the native runtime Development Cases, the current runner records the OneFLOW process exit status. The STL read report contains triangle/node/cell counts, but the GUI runner does not yet compare those fields against per-Case expected-result assertions. Therefore a green Case run currently proves successful execution of the declared runtime path; the documented node/cell counts are expected outputs to inspect, not yet an independently enforced numerical oracle. Adding machine-checked result manifests is a separate next step.

A passing GoogleTest establishes the named parser/writer behavior only.
A passing native Development Case additionally exercises runtime parameter
ingestion, Grid dispatch, topology staging and CGNS output. Neither alone
proves that an imported surface is watertight, manifold, consistently
oriented, correctly scaled, physically labelled, or suitable for CFD volume
generation. Those require independent geometry, topology, unit/provenance
and meshing contracts.

Relationship to existing Python tools
-------------------------------------

The existing grid/stl Python utilities remain useful for planar boundary
extraction and Gmsh geometry generation. They have a different responsibility:

* Python utility: planar boundary extraction, sidecar-driven domain setup,
  Gmsh generation, and independent reference checks for selected format cases.
* C++ Geometry adapter: native STL file I/O and source-facet preservation,
  integrated into the standard OneFLOW runtime and mesh conversion pipeline.
* Later topology service: explicit vertex welding, connected components,
  edge incidence, manifold checks and surface orientation diagnostics.
* Later conversion adapter: validated surface topology into a mesher or
  OneFLOW Mesh Core representation.

The Python STL development cases under ``grid/stl/python/`` currently registered in the development
case catalog exercise the Python utility only. They must not be described as
native C++ STL reader tests. Python reference cases remain useful when they
exercise distinct behavior or provide an independent cross-check; they are
not duplicate production backends. In particular, the original
``grid/stl/python/ascii_basic`` Python case is retained and registered as
``stl_python_ascii_basic`` alongside the native ``stl_ascii_basic`` Case.
They share a simple input shape but verify different implementation paths:
Python checks the reference utility's parsing and planar-boundary extraction,
while C++ checks OneFLOW's native reader, topology conversion and CGNS output. The native ``stl_ascii_basic`` and
``stl_ascii_crlf`` Cases exercise the production C++ route, with the latter
preserving CRLF line endings in its fixture. The first native C++ GoogleTest target is now
registered with CTest through the geometry test helper. Additional STL cases should be added there as reader/writer coverage expands.
The registered Python ``stl_python_ascii_basic`` and
``stl_python_ascii_crlf`` reference cases and native
``stl_ascii_crlf`` integration Case have distinct responsibilities:
the former checks the Python utility, while the latter exercises the C++
reader through OneFLOW's normal Grid task and CGNS conversion path.

Mesh Contract Workbench entry
-----------------------------

The machine-readable contract is
tools/verification/contracts/mesh/stl-native-io.json, with ID
mesh/stl-native-io. It is grouped under Mesh / Interchange Formats and
currently has status partially_verified. The Overview, Principle, Reference and
Architecture tabs describe the contract. The Implementation tab links the C++ header/source and the Verification
tab links native GoogleTests. Linked tests are declarations, not PASS results;
CI must provide execution evidence before the contract can be promoted to
verified. This prevents a source declaration from being mistaken for a
passing test.

References
----------

* Library of Congress, STL (STereoLithography) File Format:
  https://www.loc.gov/preservation/digital/formats/fdd/fdd000504.shtml
* Gmsh Reference Manual, file formats and geometry import:
  https://gmsh.info/doc/texinfo/gmsh.html
* OneFLOW mesh_architecture, gmsh_interface and stl_to_gmsh_2d pages.

Native OneFLOW grid Case
========================

``grid/stl/ascii_basic`` is a normal OneFLOW grid Case. It is registered in
``grid/case_catalog.json`` without a Python execution backend and is launched
through the standard ``OneFLOW.exe r <case-directory>`` entry point. It does not require a ``script/control.txt`` manifest. The runtime discovers
parameter files in ``script/`` according to its filename-order loading rule;
``grid.txt`` selects ``sourceGridType = "stl"`` and
``targetGridType = "cgns"``.

The C++ grid-conversion dispatcher reads the STL through ``StlMeshReader``,
constructs a two-dimensional unstructured triangle surface, and writes CGNS
through ``UnstructuredCgnsZoneWriter``. Shared coordinates are welded using
the explicit ``stlWeldTolerance`` Case parameter. The bundled two-triangle
square is expected to produce four nodes and two TRI_3 cells. A JSON read
report is written under the Case ``results`` directory. The companion
``grid/stl/ascii_crlf`` Case uses the same geometry encoded with CRLF line
endings to verify the native text parser through the production entry point.
Both native Cases omit ``control.txt`` and rely on automatic parameter file
discovery. This is a surface
conversion milestone, not a volume mesher; STL boundary-condition semantics
and physical units are not inferred.


Native OneFLOW grid Case configuration
--------------------------------------

The native `grid/stl/ascii_basic` Case is a two-dimensional surface conversion
fixture and intentionally has no `script/control.txt`. The runtime discovers
the parameter files in `script/` using its filename-order loading rule:
`dimension.txt`, `grid.txt`, and `simu.txt`. The dimension declaration is
explicit (`int dimension = 2;`); it is a runtime parameter, not an inferred
property of the STL file. The Case catalog and preflight validator also support
native Cases without a runtime manifest and validate their discovered
`script/*.txt` parameter files.


STL case directory convention
------------------------------

STL Case names describe the input and expected behavior, not the implementation
language. Native OneFLOW Cases live directly under ``grid/stl/<case>`` and use
IDs such as ``stl_ascii_crlf``. Python reference Cases live under
``grid/stl/python/<case>`` and use the ``stl_python_<case>`` catalog namespace
because catalog IDs are globally unique. This keeps the case identity readable
while making the implementation boundary explicit in the directory and
execution metadata. New STL Cases must follow this convention.


Tolerance-based vertex identity and spatial data structures
-----------------------------------------------------------

The previous native conversion and Python planar helper used coordinate-wise
rounding by the weld tolerance. Rounding is not a geometric distance test:
near points can fall on opposite cells, while points in one cell can be farther
apart than the tolerance.

The native C++ conversion now reuses the existing ADT-backed PointLocator.
The ADT is only a broad-phase spatial index; the locator applies an exact
squared Euclidean distance check in the narrow phase. The Python reference
uses a tolerance-sized spatial hash, searches the 27 neighboring buckets, and
applies the same Euclidean predicate. Both retain the first inserted point as
the representative and reuse the nearest matching representative.

This is representative-based welding, not transitive connected-component
clustering. For A-B-C where A-B and B-C are within tolerance but A-C is not,
the algorithm does not merge C merely because B matched A. This policy is
explicit because transitive closure can change topology and can chain together
points whose total spread is much larger than the tolerance.

ADT, KD-tree and uniform-grid hashing are candidate lookup structures; none
defines geometric identity. A candidate must pass the narrow-phase metric.
OneFLOW already vendors nanoflann, but a second spatial dependency is not
justified here: the existing ADT PointLocator serves the native path, while
nanoflann remains an option to benchmark if profiling later shows a need.

StlVertexWeldingTest and the Python reference tests cover cell-boundary misses,
same-cell false merges and non-transitive representative behavior. This does
not port Gmsh .geo generation, sidecar physical-boundary assignment or the
external Gmsh CLI tests into C++; those belong to a separate Gmsh conversion
adapter and must not be represented as STL parser coverage.
