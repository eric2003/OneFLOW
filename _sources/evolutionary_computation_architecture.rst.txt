Evolutionary Scientific-Computing Architecture
=================================================

Status: architectural proposal and design vocabulary. This page does not
claim that OneFLOW already has a general plugin runtime.

Objective
---------

The objective is not to build a more elaborate control-file reader. It is to
make scientific algorithms and data transformations describable, testable,
composable, replaceable, and verifiable without losing scientific meaning.

A control-file reader is a useful first specimen because it already contains
distinct representation changes:

.. code-block:: text

   legacy text
     -> syntax parsing
     -> ConfigDocument
     -> typed database entries
     -> runtime parameter state
     -> distributed state
     -> diagnostic snapshot

Each arrow needs a contract: accepted input, produced output, invariants,
side effects, failure behavior, ownership, and verification evidence.

The DNA/RNA factory analogy
---------------------------

+----------------------+--------------------------------------------------+
| Analogy              | OneFLOW meaning                                  |
+======================+==================================================+
| DNA                  | Versioned algorithm/data contracts, assumptions,|
|                      | units, and provenance                            |
+----------------------+--------------------------------------------------+
| Transcription        | Convert legacy input into a structured           |
|                      | intermediate representation                      |
+----------------------+--------------------------------------------------+
| RNA processing       | Normalize, validate, enrich, or adapt an         |
|                      | intermediate representation                      |
+----------------------+--------------------------------------------------+
| Assembly             | Compose compatible units into a pipeline or     |
|                      | solver workflow                                  |
+----------------------+--------------------------------------------------+
| Protein              | Concrete executable algorithm, data structure, |
|                      | mesh, field, or solver result                    |
+----------------------+--------------------------------------------------+
| Quality control      | Preconditions, postconditions, invariants,      |
|                      | tests, verification, and provenance checks      |
+----------------------+--------------------------------------------------+
| Mutation/evolution   | Reviewed changes to algorithms, schemas,        |
|                      | adapters, and implementations                    |
+----------------------+--------------------------------------------------+

A small algorithm robot is not necessarily a class, plugin, or dynamically
loaded library. It is first a unit of responsibility with an explicit
contract. A free function can qualify if its inputs, outputs, state effects,
and tests are clear. Do not introduce a factory or plugin manager merely to
complete the analogy.

Machine-readable computation-unit contract
-------------------------------------------

The proposed descriptor is defined in
``doc/specs/computation-unit.schema.json``. It records identity/version, role,
typed input and output ports, units, ownership, preconditions, postconditions,
side effects, failure policy, determinism, lifecycle, compatibility adapters,
verification evidence, known limits, and scientific provenance.

This schema is a design contract, not an executable runtime ABI. Its first
purpose is to make architecture review precise and machine-checkable. Later
tools may validate descriptors, build dependency graphs, generate test
matrices, or render pipeline diagrams. Tooling must not infer missing
scientific assumptions.

Intermediate representations
-----------------------------

Avoid direct coupling between every input format and every consumer. Prefer
a small number of explicit intermediate representations (IRs), each with a
documented semantic scope. Candidate representations include:

* ConfigDocument: parsed configuration syntax before runtime commit.
* Typed parameter snapshot: validated values with declared types/defaults.
* Mesh snapshot: canonical geometry/topology values independent of solver globals.
* Topology and boundary views: derived connectivity with invalidation rules.
* Algorithm/workflow description: selected algorithms, dependencies, and scope.
* Result/diagnostic snapshot: outputs with units, provenance, and verification.

An IR is useful when multiple producers or consumers need to interoperate,
or when it freezes a meaningful semantic boundary. It is not useful if it
merely duplicates a data structure. Adapters should be small and directional.
A legacy input adapter should not own solver lifecycle or numerical policy.

Architecture principles
-----------------------

Representation and meaning
~~~~~~~~~~~~~~~~~~~~~~~~~~

Parsing answers what was written. Validation answers whether it is well-formed
and type-correct. Domain interpretation answers what it means for this case.
Application answers which runtime state changes. Execution answers what
computation is performed. Keep these contracts distinguishable.

Explicit state and side effects
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A component that reads a file, mutates a global database, performs a parallel
collective, or changes device state is not a pure transformation. Declare
the effect and lifecycle scope. Prefer value snapshots at boundaries. Isolate
legacy globals behind adapters and explicit setup/teardown contracts.

Replaceable implementations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Keep stable external contracts while moving implementations behind module
boundaries. Version schemas and IRs when meaning changes. Use compatibility
adapters for old representations instead of scattering version checks
throughout numerical algorithms.

Composable failure behavior
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Every stage must define whether it leaves state unchanged, rolls back its own
mutations, or leaves diagnostics. Cleanup ownership must be unambiguous.
Case-scoped operations must not finalize process-scoped resources. A failed
stage must not silently contaminate the next case.

Tests as executable architecture constraints
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Unit tests pin local invariants.
* Component tests pin collaboration between adjacent transformations.
* Integration tests pin the real production pipeline and its side effects.
* Regression tests pin scientific behavior against known cases.
* Manufactured-solution, conservation, and convergence tests provide evidence
  about numerical correctness.

If a small contract test requires linking the entire solver executable,
inspect the dependency boundary. This is not automatic proof that another
abstraction is required. A green test establishes only the properties it asserts.

Scientific provenance and verification
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: text

   published theory/reference
     -> OneFLOW mathematical formulation
     -> implementation and configuration
     -> focused property tests
     -> manufactured/benchmark case
     -> regression evidence and known limits

Citations establish provenance, not correctness. Unit tests establish local
properties, not universal model validity. Numerical baselines must be
reproducible and reviewed before being changed.

Current specimen: control-file ingestion
-----------------------------------------

.. code-block:: text

   script/control.txt
     -> ControlFileReader: resolve ordered script list
     -> ConfigLoader: parse legacy text
     -> ConfigDocument: structured representation
     -> ConfigDatabaseAdapter: commit typed entries
     -> DataPara: runtime parameter storage
     -> ControlDatabaseIO: broadcast and diagnostic dump

The source responsibilities are now separated across five production units:

* ``ControlFileReader.cpp`` coordinates manifest reading and ordered loading.
* ``ControlScriptParser.cpp`` adapts file/parser input to ConfigLoader and
  commits the resulting ConfigDocument.
* ``ControlParameterResolver.cpp`` resolves literal or database-backed array
  lengths required by legacy syntax.
* ``ControlDatabaseIO.cpp`` owns database broadcast, compression, and dumps.
* ``ParaFile.cpp`` retains older parameter inspection helpers.

The public declarations in ``ParaFile.h`` remain stable. The parser adapter
uses the existing runtime resolver through a narrow translation unit rather
than through the broad legacy inspection module. The existing ConfigLoader
unit tests exercise injected resolver callbacks; a new resolver test target
exercises literal lengths and named lookup through the real DataBase. This is
still not proof of a complete multi-file reader integration path.

Expected problems and response strategy
----------------------------------------

+------------------------------+------------------------------------------+
| Problem                      | Response                                 |
+==============================+==========================================+
| New source not compiled      | Inspect CMake discovery and build graph  |
+------------------------------+------------------------------------------+
| Reader still links legacy    | Extract parser adapter after dependency  |
| parser                       | mapping; never duplicate parser logic    |
+------------------------------+------------------------------------------+
| Second case inherits state   | Run A then B in one process; assert B's  |
|                              | own defaults and teardown contract       |
+------------------------------+------------------------------------------+
| Reader hangs or reads none   | Initialize the real serial backend       |
+------------------------------+------------------------------------------+
| Cleanup hides original error | Preserve first exception and assert      |
|                              | cleanup                                   |
+------------------------------+------------------------------------------+
| Plugin framework grows early | Start with ordinary functions and values |
+------------------------------+------------------------------------------+
| Duplicate IRs proliferate    | Require real producer/consumer decoupling |
+------------------------------+------------------------------------------+
| Metadata drifts              | Add schema validation and traceability   |
+------------------------------+------------------------------------------+
| Unit tests mistaken for      | Separate software tests from numerical   |
| scientific validation        | evidence                                  |
+------------------------------+------------------------------------------+

Evolution roadmap
-----------------

Short term: make today's boundaries real
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#. Complete the control-reader split without changing parameter semantics.
#. Confirm automatic source discovery and compile on Windows and Linux.
#. Extract only parser adapter dependencies needed to link the real reader.
#. Build a temporary multi-file fixture; initialize the real serial runtime;
   assert typed values and log/database.log.
#. Add sequential-case and failure-cleanup tests at the correct lifecycle boundary.
#. Keep the handoff, architecture page, and integration plan aligned with evidence.

Medium term: establish reusable data contracts
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#. Stabilize ConfigDocument and typed parameter snapshots as versioned IRs.
#. Apply the same pattern to mesh snapshots, topology, boundary data, and results.
#. Define ownership, units, invalidation, and lifecycle rules.
#. Add schema validation, dependency diagrams, and test-to-contract traceability.
#. Introduce compatibility adapters where legacy APIs must coexist.

Long term: an algorithm and dataflow ecosystem
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#. Describe solver pipelines as compositions of compatible computation units.
#. Allow alternative implementations to satisfy the same mathematical contracts.
#. Support CPU, accelerator, serial, and distributed implementations through
   declared capabilities and effects.
#. Generate test plans and provenance reports from descriptors and pipeline graphs.
#. Compare algorithm variants using reproducible verification suites.
#. Consider dynamic plugin loading only when real extension requirements justify it.

The target is not a universal abstraction layer over every line of code. It is
a small set of stable semantic contracts that lets the system change shape
without losing correctness, traceability, or scientific meaning.

Review questions for each future refactor
-----------------------------------------

#. Which producer and consumer are being decoupled?
#. What representation crosses the boundary, and what does it mean?
#. Which invariants, units, ownership, and lifecycle rules apply?
#. What side effects occur, and who owns rollback?
#. Which test proves the contract, and what higher-level evidence is missing?
#. Can an alternative implementation be substituted without changing consumers?
#. Is this a reusable boundary, or merely another wrapper?

If these questions have no concrete answers, keep the implementation simple
and defer the abstraction.


The resolver's focused tests are registered in `tests/io/CMakeLists.txt` as `control_parameter_resolver_test`. They intentionally link the real DataBase implementation instead of mocking named-parameter lookup. This pins the adapter contract without claiming that the full control manifest, broadcast, and diagnostic pipeline has been integrated.
