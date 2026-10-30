OneFLOW Architecture Priority Roadmap
=======================================

Purpose
-------

This roadmap ranks work by system-wide return, not by the number of files
changed or by the apparent size of a subsystem. It is intended to prevent
local refactoring from becoming an end in itself.

The ranking is a working order, not a claim that lower-ranked areas are
unimportant. Reassess it when production evidence changes the risk profile.

How to rank work
----------------

Use this qualitative model when comparing candidate changes:

.. math::

   P = \frac{I \cdot F \cdot B \cdot U}{C + V}

where:

* :math:`I` is the impact if the capability is correct;
* :math:`F` is how frequently the capability is exercised;
* :math:`B` is the blast radius when it fails;
* :math:`U` is uncertainty or the lack of trustworthy evidence;
* :math:`C` is implementation cost;
* :math:`V` is the cost of building convincing verification evidence.

This is a prioritization aid, not a measured metric. Do not assign false
precision to the factors. In particular, high uncertainty increases the
value of an investigation or a focused integration test; it does not
automatically justify a large rewrite.

Ranked investment sequence
--------------------------

.. list-table::
   :header-rows: 1
   :widths: 7 27 23 28 15

   * - Rank
     - Investment
     - System-wide leverage
     - Required evidence
     - Current posture
   * - 1
     - Runtime lifecycle, configuration ingestion, and case isolation
     - Every solver run depends on correct setup, state ownership, cleanup,
       and reproducible case paths.
     - Real entry-point tests; multi-file configuration; sequential cases;
       setup and execution failure cleanup; Windows/Linux CI.
     - Active: R1 partial.
   * - 2
     - Mesh core, topology, boundary, and derived-data lifecycle
     - Geometry, conversion, partitioning, solver connectivity, visualization,
       and mesh workflows share these semantics.
     - Invariant tests; identity/orientation tests; cache invalidation;
       representative 1D/2D/3D fixtures; independent snapshots.
     - Next major engineering mainline.
   * - 3
     - Numerical algorithms and scientific verification
     - Numerical correctness determines whether successful execution
       produces scientifically meaningful results.
     - Derivation and assumptions; manufactured solutions; conservation;
       convergence/order studies; reference comparisons; real-case regression.
     - Build a traceable algorithm inventory and deepen one method at a time.
   * - 4
     - Computation-unit contracts and dataflow composition
     - Enables algorithms and data representations to be replaced or composed
       without coupling every producer to every consumer.
     - Typed ports, units, ownership, pre/postconditions, side effects,
       failure policy, provenance, compatibility tests.
     - Establish contracts and adapters before runtime machinery.
   * - 5
     - Workflow and GUI developer experience
     - Improves usability, repeatability, artifact dependencies, and diagnosis.
     - End-to-end workflow tests; artifact provenance; restart/stop semantics;
       no source-tree pollution.
     - Advance in step with stable runtime and mesh contracts.
   * - 6
     - Documentation and developer-tool automation
     - Reduces maintenance cost and makes architecture discoverable across
       all other investment areas.
     - Documentation build; API extraction; link and formula checks;
       diagrams kept consistent with implementation.
     - Continuous cross-cutting work, not a substitute for implementation.

Why this order
--------------

Runtime correctness is the first gate: an algorithm cannot be trusted if
the case parameters, global state, paths, or device/parallel resources are
not initialized and released consistently.

Mesh is the next broad leverage point. A mesh is not merely an input file or
a collection of coordinates. It carries geometry, entity identity, oriented
connectivity, boundary semantics, and derived caches. A stable mesh core
therefore improves many independent workflows at once.

Numerical methods follow as a distinct scientific-quality investment. The
implementation must be connected to its mathematical model and verification
evidence. A passing software test suite does not by itself establish
consistency, stability, conservation, or convergence of a numerical method.

Computation-unit contracts are an enabling architecture. They should emerge
from stable boundaries found in runtime, mesh, and numerical work. Do not
build a generic plugin host before real algorithms need independently
replaceable implementations.

Workflows and GUI then consume those stable capabilities. Documentation
automation supports every stage, but the authoritative design narrative,
scientific derivations, and system diagrams remain curated documentation.

Current mainline: R1 runtime integration
----------------------------------------

The current configuration path can be described as a sequence of semantic
transformations and effects:

.. code-block:: text

   case directory
        |
        v
   script/control.txt
        |
        v
   ordered script files
        |
        v
   legacy syntax parser
        |
        v
   ConfigDocument (parsed representation)
        |
        v
   typed DataBase commit
        |
        v
   parallel distribution / diagnostic snapshot
        |
        v
   CaseIOSettings validation
        |
        v
   project and domain I/O roots
        |
        v
   io_paths.txt and case execution

The arrows are not interchangeable. Parsing determines what the text says;
typed commit determines how values enter the runtime database; validation
determines whether case-specific settings are acceptable; application mutates
project path state; parallel distribution and diagnostics have process and
filesystem effects.

The next R1 acceptance sequence is:

#. Exercise the production `ReadControlInfo(caseDir)` with a temporary
   multi-file manifest and a correctly initialized serial parallel runtime.
#. Assert typed values in the real database and verify the database diagnostic
   output.
#. Exercise `SetupCaseEnvironment()` far enough to verify effective project,
   mesh, and flow roots and the generated path diagnostics.
#. Run two distinct cases sequentially and prove that the second case cannot
   observe the first case's case-scoped settings.
#. Inject a setup failure and verify case cleanup while preserving the
   original execution exception.
#. Keep process-scoped finalization owned by the process lifecycle; do not
   finalize shared parallel/device resources from case teardown.

A focused parser test, resolver test, or path helper test proves only its
declared contract. It must not be reported as proof of this whole sequence.
See :doc:`runtime_architecture` and
`doc/plans/setup-case-environment-integration-test.md` for source-level
details and open questions.

Mesh mainline after R1
---------------------

The mesh investment should preserve a strict distinction between source data,
canonical semantic snapshots, and derived views:

.. code-block:: text

   mesh reader / generator
           |
           v
   canonical mesh snapshot
      /      |       \
     v       v        v
   geometry topology boundary
              |        |
              v        v
          derived connectivity
                 |
                 v
          solver / partition /
          visualization adapters

The canonical snapshot should not depend on solver globals or GUI state.
Topology and boundary views must document orientation, identity, ghost-entity
rules, and cache invalidation. Serialization (JSON, VTU, or another format)
should consume snapshots; it should not reach back into live solver state to
reconstruct meaning.

Numerical algorithm work
------------------------

For each algorithm selected for deeper work, maintain a traceable chain:

.. code-block:: text

   physical model and assumptions
              |
              v
   mathematical formulation
              |
              v
   discrete equations and boundary treatment
              |
              v
   implementation and data contract
              |
              v
   focused unit / component tests
              |
              v
   manufactured-solution and conservation checks
              |
              v
   grid convergence and reference-case comparison

For a grid-refinement study with characteristic spacing :math:`h`, an
observed order can be estimated from errors :math:`E_h` and
:math:`E_{h/r}`:

.. math::

   p_{obs} = \frac{\ln(E_h/E_{h/r})}{\ln(r)}.

This estimate is meaningful only when the solutions are in the asymptotic
convergence range, the error norm and refinement ratio are defined, and the
compared runs use consistent physical and numerical settings. The document
for a specific method must state these assumptions rather than presenting
the formula as a universal pass/fail test.

Documentation and source comments
----------------------------------

Use two complementary documentation layers.

**Source-level API documentation** should be generated from useful, stable
comments (for example, Doxygen-compatible comments) on public classes,
functions, data structures, and non-obvious invariants. Comments should state
intent, units, ownership/lifetime, preconditions, postconditions, side effects,
and failure behavior where these are part of the contract. Avoid narrating
obvious statements or copying entire design documents into headers.

**Curated documentation under `doc/source`** remains necessary for material
that spans files or requires explanation: architecture decisions, subsystem
maps, mathematical derivations, discretization choices, algorithm limits,
workflow diagrams, numerical experiments, reference comparisons, and
migration plans. Generated API pages should be linked from these guides; they
should not replace them.

Use reStructuredText/Sphinx capabilities consistently:

* `.. math::` for equations and derivations;
* `.. figure::` for explanatory diagrams and annotated plots;
* `.. graphviz::` or a maintained text diagram when the build configuration
  supports it;
* `.. list-table::` for contracts, test matrices, and comparison summaries;
* citations and a reference list for external scientific claims.

A diagram must explain a real boundary, data transformation, lifecycle, or
algorithm. Keep it adjacent to the prose that interprets it. Every diagram
should have a caption and should not imply that a planned component already
exists. If an automatically generated diagram is sourced from code, identify
the generator and keep a reviewed conceptual diagram for the system-level
story.

Review checkpoint
-----------------

At each substantial milestone, answer:

* Which high-leverage system risk has been reduced?
* Which users, algorithms, or subsystems benefit from the change?
* What behavior is now verified, and by which exact test or CI revision?
* Which important behavior remains unverified?
* Did the evidence change the priority ranking?
* Did the change introduce an abstraction that has no second real use case?

The goal is not maximum abstraction or maximum test count. The goal is a
scientific-computing system whose state, algorithms, data transformations,
and numerical claims can be understood, composed, replaced, and verified.
