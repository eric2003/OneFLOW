Contract Workbench Information Architecture
============================================

Purpose
-------

The Contract Workbench is a shared information architecture for engineering
knowledge and verification evidence. It is not a formula workbench with other
domains forced into formula-shaped fields. Formula Knowledge is one domain
among several; Case, Runtime, Workflow, Mesh, Geometry and Algorithm have
different native concepts and evidence.

The common structure provides predictable navigation and traceability while
each domain retains a typed payload and domain-specific visualizers.

Shared navigation
-----------------

.. list-table::
   :header-rows: 1
   :widths: 20 42 38

   * - Section
     - Question answered
     - Typical content
   * - Overview
     - What is it, what problem does it solve, and what is its current status?
     - Intent, scope, owner, domain/kind, status, claims and limitations.
   * - Principle
     - Why is it designed this way?
     - Physical or mathematical theory, engineering rationale, invariants,
       assumptions, lifecycle rules and limiting cases.
   * - Reference
     - Which external or internal authorities inform it?
     - Papers, standards, software manuals, upstream implementations, ADRs
       and related contracts.
   * - Architecture
     - Where does it sit in OneFLOW and what interacts with it?
     - Module boundaries, dependencies, ownership, data flow, state diagrams
       and sequence diagrams.
   * - Implementation
     - How is the behavior realized?
     - Algorithms, data structures, lifecycle, configuration, adapters,
       execution steps and important invariants.
   * - Code
     - Where is the implementation in the repository?
     - Source paths, symbols, call chains, bounded source excerpts and tests.
       Code excerpts must be traceable to repository files.
   * - Documentation
     - Where is the maintained explanation or user guidance?
     - Sphinx pages, developer guides, user manuals and design documents.
   * - Verification / Evidence
     - What has actually been checked, and what remains unverified?
     - Test targets, exact commands/configuration, case runs, artifacts,
       platform/revision, outcome, freshness and explicit evidence limits.

These sections are common navigation, not a requirement that every contract
must contain every kind of content. Empty or undeclared information must be
shown as such; the UI must not invent a formula, diagram, reference or PASS
state merely to fill a tab.

Domain-specific knowledge and visualizers
-----------------------------------------

Shared navigation does not erase domain distinctions. Each domain has a
typed payload and may add focused views:

.. list-table::
   :header-rows: 1
   :widths: 18 42 40

   * - Domain
     - Native knowledge structure
     - Optional domain-specific views
   * - Formula Knowledge
     - Equations, symbols, notation aliases, assumptions, dimensional
       consistency, derivations and limiting behavior.
     - Formula / LaTeX preview, theory curves, worked examples.
   * - Parameter
     - Parameter identity, aliases, type, units, source/load order, default,
       override behavior, consumers and runtime read timing.
     - Formula association, runtime consumers, examples, declaration history.
   * - Case
     - Physical setup, solver/model, mesh, numerical methods, initial and
       boundary conditions, expected outputs and comparison policy.
     - Case setup, input files, suite membership, run artifacts and baselines.
   * - Runtime
     - Process/task state machine, resource ownership, transitions, failure
       handling, cancellation, recovery and concurrency constraints.
     - State-transition diagram, event timeline, resource/lifecycle trace.
   * - Workflow
     - Steps, dependencies, artifact contracts, input/output bindings,
       scheduling and failure propagation.
     - Workflow graph, run history, artifact lineage and diagnostics.
   * - Mesh
     - Format-neutral entities and topology, coordinates, connectivity,
       boundaries, regions, identity, orientation and quality; plus Gmsh,
       Plot3D, CGNS and future adapter semantics.
     - Grouped Mesh Knowledge catalog, topology inspector, format mapping,
       boundary map, quality metrics and round-trip comparison.
   * - Geometry
     - Geometric entities, CAD topology, construction history and geometry-to-mesh
       association.
     - Geometry scenarios, entity graph, topology fixtures and geometry checks.
   * - Algorithm
     - Mathematical model, discretization, data dependencies, algorithm
       stages, complexity, stability/conservation properties and assumptions.
     - Pseudocode, dependency graph, convergence/order plots and benchmarks.

Formula, curve, mesh, workflow graph and case setup are domain visualizers,
not mandatory shared tabs. A mesh contract should show topology and boundary
semantics rather than a fake equation; a Runtime contract should show a state
machine rather than a curve; an Algorithm contract may show both equations
and a data-flow diagram where each is relevant.

Shared contract envelope and typed payload
-------------------------------------------

The shared envelope should remain small and versioned:

* stable contract ID and schema version;
* domain, kind, title, status and tags;
* overview, claims and limitations;
* links to principle, references, architecture, implementation and
  documentation;
* typed evidence declarations and their resolution/execution state;
* a domain-specific payload validated by its own schema.

Domain details belong in the typed payload rather than an ever-growing set of
unrelated optional fields. The existing Formula Knowledge catalog and
geometry scenario fixtures remain authoritative in their own domains while
adapters expose their shared navigation metadata. Do not duplicate their
editable source data into a second catalog during migration.

Evidence semantics
------------------

Evidence is a bounded claim, not a decoration. Each evidence record should
identify its kind, target, exact command and configuration where applicable,
source revision, platform/backend, outcome, artifact location and limitations.
Useful evidence kinds include schema/static checks, unit/analytical tests,
subsystem integration, production-path integration, representative case
regressions, cross-platform/backend checks and independent benchmarks.

Keep declaration and execution separate. A linked test path means the test is
declared; it does not mean it ran. A passing unit test does not establish
physical accuracy, and a successful mesh parse does not establish semantic
boundary preservation. Display not-run, stale, blocked, failed and passed
states explicitly, with evidence scope and freshness.

Architecture and ownership
---------------------------

.. code-block:: text

   Domain sources / typed payloads
      Formula catalog | Case contracts | Runtime contracts
      Workflow definitions | Mesh scenarios | Algorithm records
                         |
                         v
             Shared contract envelope
                         |
                         v
       Qt-independent catalog / validation
                         |
             +-----------+-----------+
             |                       |
      Evidence resolver        Typed domain adapters
             |                 (case, runtime, mesh, ...)
             +-----------+-----------+
                         |
                         v
             Shared Workbench navigation
                         |
                         v
        Domain-specific views and evidence

Catalog validation, evidence resolution and domain logic must remain Qt
independent. GUI code presents data and invokes explicit adapters; it must
not duplicate solver algorithms or infer verification outcomes. Execution
continues to use the authoritative CTest, Python unittest, CFD/grid, Case and
Workflow services.

Current implementation status
------------------------------

The GUI Contract Workbench currently exposes the common Overview, Principle,
Reference, Architecture, Implementation, Code, Documentation and
Verification / Evidence navigation for the typed contract domains. Formula
and Parameter retain formula and parameter-specific views; Case and Workflow
retain their setup/graph views; Geometry continues to use its fixture-backed
topology explorer. The common tab layout is a navigation milestone, not a
claim that every domain's content has been fully populated.

Remaining work is to improve domain payload coverage, especially Mesh and
Algorithm; connect production-backed mesh evidence; enrich Runtime and
Workflow execution evidence; and audit each tab for source-backed content.
Each new claim must distinguish implemented, declared, tested and not yet
verified behavior.

Related documents
-----------------

* :doc:`workbench_architecture`
* :doc:`workbench_user_guide_zh`
* :doc:`mesh_architecture`
* :doc:`gmsh_interface`
* :doc:`case_parameter_design`
