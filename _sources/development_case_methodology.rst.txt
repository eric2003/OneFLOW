Development Case Methodology
=============================

Purpose
-------

OneFLOW Development Cases are reusable engineering examples, not merely
small regression tests. They demonstrate how to formulate a mesh or CFD
problem, construct its inputs, invoke the relevant implementation, transfer
artifacts between stages, and establish evidence that the result satisfies a
specific contract. A well-designed Case should teach a future contributor
how to build the next Case without copying accidental details from an
unrelated example.

This guide defines the common method for 2D and 3D mesh-development Cases.
It complements the executable catalog in ``grid/case_catalog.json``, the
source tests, the Workbench selectors, and the workflow evidence records.

Case design principles
----------------------

Each Case must answer one bounded engineering question. Its name and
acceptance criteria must say exactly what it demonstrates. Do not label a
Case "two connected blocks" if the shared interface has extra interior nodes,
or "automatic triangles" if the element count is hard-coded.

Use a progressive sequence:

* establish the smallest valid geometry and topology;
* add one new topological or algorithmic property;
* add automatic or more general construction only after the simpler
  contract is independently verified;
* add multi-block, curved, internal-boundary, or 3D complexity as separate
  Cases rather than hiding it in a minimal example.

A Case is complete only when its inputs, procedure, outputs, limitations,
tests, documentation, and GUI discovery agree.

2D and 3D organization
---------------------

The Case Catalog carries optional ``dimension`` metadata with value 2 or 3.
The GUI presents a dimension-first hierarchy, followed by the Case family and
individual Case. The native C++ and Python reference variants are siblings
under the same family; they are not separate conceptual categories.

Recommended source organization:

.. code-block:: text

   grid/
     mesh-layout/
       2d/<case_name>/
       3d/<case_name>/
     stl_to_gmsh/
       <geometry_or_problem>/
         cpp/
         python/

The existing repository contains legacy paths and domain-specific families;
do not relocate a Case solely to make the tree look uniform. New Cases should
follow the dimension-first convention where it improves discoverability.
A missing or ambiguous dimension belongs in the Unclassified group until
the actual mesh/layout/profile establishes it. Do not infer dimension from
the title alone.

A 2D surface mesh and a 3D volume mesh are different contracts. A 3D STL is
usually a surface representation used to construct a volume mesh; its
triangle facets do not make the resulting volume mesh 2D. Record both the
input representation and the target mesh dimension where that distinction
matters.

Recommended Case contents
-------------------------

A producer Case should keep only the files needed to describe and reproduce
that Case:

.. code-block:: text

   <case>/
     README.md
     input/                 # immutable source geometry/data
     script/                # runtime parameters and construction layout
     cpp/                   # only when source/build glue is Case-specific
     python/                # only for an independent reference implementation
     expected.json          # Python reference output contract, when applicable

Use the actual existing family conventions when they differ. Do not add empty
folders or duplicate shared scripts just to satisfy a template. Shared
algorithms belong in shared implementation modules; per-Case folders should
contain declarative inputs and only genuinely Case-specific code.

Construction and execution sequence
-----------------------------------

Document the following sequence for each mesh producer:

#. **Question and scope** — what single capability is being demonstrated,
   and what is deliberately outside its scope?
#. **Input model** — geometry source, units, coordinate frame, topology,
   boundary labels, and whether the input is shared with a reference Case.
#. **Construction** — how geometry entities, curves, surfaces, blocks,
   transfinite constraints, mesh sizes, and physical groups are created.
#. **Generation** — the exact entry point, configuration files, external
   tools, versions/options that affect reproducibility, and the expected
   artifact path.
#. **Inspection** — nodes, element types/counts, boundary elements,
   connected components, manifold/incidence checks, orientation, and
   geometric quality as relevant.
#. **Conversion** — how the artifact is handed to CGNS and OneFLOW, including
   boundary-condition mapping and identity preservation.
#. **CFD use** — whether a solver is actually run, and whether that run is a
   smoke test, regression, verification, or accuracy benchmark.
#. **Evidence** — link every claim to the relevant unit/integration test,
   workflow result, artifact, or numerical reference.

Do not claim every Case must run a CFD solve. A single-cell topology fixture
may correctly stop after mesh import and boundary-map verification. Conversely,
a Case intended to establish numerical accuracy needs more than a successful
mesh conversion or a two-step solver smoke test.

Native C++ and Python reference responsibilities
-------------------------------------------------

When both implementations are provided:

* use the same immutable geometry input and equivalent declarative settings;
* keep generated artifacts and reports outside the source tree;
* define the artifact basename and semantic boundary names consistently;
* compare invariants and acceptance criteria, not necessarily node numbering,
  element ordering, or unstructured connectivity;
* keep Python implementation independent enough to serve as a useful
  reference, rather than calling the C++ implementation under test;
* record implementation-specific limitations and external dependencies.

The Python ``expected.json`` should express stable properties: element
families, minimum counts, named Physical Groups, bounds, quality thresholds,
or expected artifact names. Avoid exact counts when the mesher is explicitly
in automatic mode and the count is allowed to vary with mesh size or Gmsh
version.

Test and validation layers
--------------------------

Separate the evidence by what it proves:

* **Unit tests** — parsing, layout interpretation, topology checks, naming,
  parameter validation, and error diagnostics.
* **Producer integration tests** — run the real C++ or Python producer and
  inspect the resulting mesh artifact.
* **Conversion integration tests** — Gmsh to CGNS and CGNS to OneFLOW,
  including element types, boundary locations, BC types, region maps, and
  grid metadata.
* **Mesh-to-CFD workflow tests** — prove that the producer's artifacts are
  correctly bound into the next stages. A solver smoke test proves runtime
  compatibility only, not numerical accuracy.
* **Numerical verification** — compare against analytic solutions,
  manufactured solutions, trusted references, convergence rates, and
  conservation/error metrics when the Case claims numerical correctness.

Tests belong in the test hierarchy (for example ``tests/python/grid`` and
the configured C++ test targets), not beside immutable source geometry when
they create runtime outputs. Every executable Case must be staged under the
user/build workspace. Never run a checked-in Case in-place if its relative
paths can create files in the repository.

Catalog and GUI registration
----------------------------

Register each runnable implementation in ``grid/case_catalog.json`` with
a unique ID, descriptive title, source directory, family ``group``, and
spatial ``dimension`` when known. A Python reference entry also declares
its execution kind and shared runner script. Keep the conceptual Case name
and implementation variant visible in the title, but avoid encoding every
configuration option into the name.

The GUI is a selector, not the source of truth. It should read the catalog
and show:

* dimension: 2D, 3D, or Unclassified;
* family: for example STL, STL to Gmsh, or Mesh Layout;
* Case: the actual runnable producer or fixture.

Mesh-to-CFD Validation has a narrower list of supported mesh producers. It
should reuse the same dimension metadata but only expose Cases for which a
valid artifact-bound conversion graph exists. Do not list arbitrary
development tests there simply because they produce a file.

Documentation requirements
--------------------------

Each Case README should explain its own inputs, construction choices,
execution command, output contract, acceptance criteria, limitations, and
where its automated tests live. Family-level documents should explain the
progression and shared conventions. This methodology document should explain
the reusable engineering process, not repeat every Case's data.

Update documentation in the same functional change as the implementation.
Keep claims synchronized across Case README, family README, Sphinx docs,
catalog metadata, tests, and verification contracts. Mark unverified claims
as such rather than implying they passed.

Review checklist
----------------

Before merging a new Case into the main development branch, check:

* [ ] one clearly bounded question and honest name;
* [ ] dimension and family metadata are correct;
* [ ] inputs are minimal, immutable, and documented;
* [ ] C++ and Python variants share the same source data where applicable;
* [ ] output paths are isolated from the repository;
* [ ] artifact names and boundary semantics are stable;
* [ ] tests cover both expected behavior and invalid inputs;
* [ ] downstream conversion claims are supported by integration evidence;
* [ ] GUI registration and selector grouping are correct;
* [ ] documentation and verification-contract links are complete;
* [ ] limitations and unverified capabilities are explicit.


3D cube progression
-------------------

The STL-to-Gmsh cube family is organized under
``grid/stl_to_gmsh/cube/<strategy>/{cpp,python}``. The three strategies
share one byte-identical cube STL: ``single_hex`` (one HEXA_8),
``six_tetrahedra`` (fixed six TETRA_4), and
``automatic_tetrahedra`` (Gmsh-controlled tetrahedral refinement).
This sequence separates a fixed topology contract from automatic meshing;
the latter must not assert a fixed element count.
