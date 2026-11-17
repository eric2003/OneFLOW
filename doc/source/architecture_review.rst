OneFLOW Architecture Review and Development Principles
======================================================

Purpose and status
------------------

This chapter is the project-wide architecture checkpoint for developers and AI
assistants. It complements the runtime, mesh, and Workbench design pages. Those pages
explain a subsystem in detail; this page asks whether the overall system is progressing
toward a coherent engineering application, what major gaps remain, and whether the next
local change addresses a meaningful part of that picture.

It separates three kinds of knowledge:

* **General CFD/CAE concerns** are a checklist of areas that mature engineering software
  commonly needs to address. Listing an area does not mean OneFLOW implements it.
* **OneFLOW observations** are source-backed findings or explicitly recorded gaps.
  Recheck code and tests when the branch changes.
* **Proposed direction** describes possible architecture, not existing capability.

This is a living review framework, not a claim that the architecture is complete.

Highest-value improvement sequence
---------------------------------

This ordering is a working risk assessment, not a claim that every item is equally
urgent or already implemented. Reorder it only when source inspection, tests, or
representative runs provide stronger evidence.

1. **Runtime lifecycle and state isolation.** Close the real case setup/read/apply/
   teardown loop, especially failure cleanup and repeated cases in one process. The
   runtime bridges explicit contexts to legacy global database, path, solver, and
   accelerator state. This boundary can invalidate otherwise correct individual
   modules and is the current active R1 slice.
2. **Mesh topology and derived-data validity.** Specify ownership and invalidation for
   node/face/cell connectivity, boundary/interface maps, ghost data, and cached
   cell-to-cell/cell-to-face relations. Stale derived data can affect many algorithms
   while appearing locally plausible. Prefer invariants and mutation-lifecycle tests
   before broad data-structure replacement.
3. **Numerical verification as a traceable chain.** For high-impact solver algorithms,
   connect governing equations and assumptions to source formulas, limiting cases,
   manufactured/analytic tests, conservation, convergence, and established benchmark
   cases. Neither unit tests nor citations alone constitute physical validation.
4. **Configuration schema and artifact provenance.** Keep parsing, type/schema
   validation, semantic resolution, and application distinct; make workflow inputs,
   outputs, staging roots, and producer-consumer dependencies explicit. This reduces
   failures caused by ambiguous paths and implicit file handoff.
5. **Build/test evidence and reproducibility.** Make focused targets, CTest discovery,
   platform differences, and exact-commit CI evidence easy to reproduce. A larger test
   count has value only when assertions cover the intended contracts.
6. **Documentation generation and interface drift.** Generate mechanical API details
   from structured source comments where stable, but keep architecture, rationale,
   equations, assumptions, diagrams, and evidence curated in ``doc/source``. Start
   with one documented module and measure generated-output stability before expanding.

This sequence prioritizes state validity and numerical trust over cosmetic renaming
or broad comment expansion. Improve names and comments as part of the relevant
boundary work, not as a competing project.

Review principles
-----------------


1. **Keep the whole-system model visible.** Before a local refactor, identify the
   subsystem, upstream inputs, downstream consumers, state ownership, lifecycle, and
   relevant user workflow.
2. **Solve a real problem, not a pattern-shaped problem.** Introduce a class, interface,
   registry, adapter, or abstraction only when it improves a concrete contract, test seam,
   dependency boundary, or maintenance cost.
3. **Move between zoom levels.** Work at the smallest useful code boundary, then return
   to the subsystem and whole-system views. A locally clean implementation can still
   preserve an undesirable system dependency.
4. **Treat tests as architecture evidence.** Focused tests state local contracts;
   integration tests show composition; numerical regressions protect engineering results.
   None substitutes for the others.
5. **Make ownership and failure paths explicit.** Track who creates, owns, mutates,
   persists, and releases important objects and artifacts, including partial setup and
   failed runs.
6. **Separate facts from plans.** Current source and tests establish behavior. Design
   pages explain intended boundaries. Roadmaps do not prove implementation.
7. **Prefer incremental compatibility.** Record behavior that must remain stable before
   moving boundaries. Avoid broad rewrites when a smaller change can test the same idea.
8. **Return to the system map at milestones.** Ask whether the overall architecture has
   improved, whether the biggest risk has changed, and whether the next task should move.

Design intent as an executable contract
---------------------------------------

AI can generate implementation quickly, but generation is only reliable when the intended
behavior and boundaries are explicit. OneFLOW should treat design intent as the source of
constraints, not as a substitute for code or tests.

Use four linked levels:

* **Intent:** explain why a capability exists, which responsibility it owns, and which
  responsibilities remain outside its boundary.
* **Contract:** define inputs, outputs, ownership, lifecycle, invariants, failure behavior,
  and compatibility requirements precisely enough to review.
* **Implementation:** write or generate code that realizes the contract. Generated code is
  still ordinary code and must remain readable, reviewable, and maintainable.
* **Evidence:** add focused tests for the contract, integration tests for composition, and
  representative numerical regressions where engineering behavior is involved.

Documentation-driven code generation is appropriate when the source description is
structured and deterministic: field tables, enum mappings, command schemas, repetitive
registrations, serialization layouts, or boilerplate derived from a stable interface.
Keep the generator and its input schema versioned, make output reproducible, and verify
that regenerating produces no unexplained diff.

Do not generate an algorithm from a broad prose statement such as "make the limiter robust."
For numerical methods, prose should first identify the governing equation, assumptions,
boundary treatment, stability or monotonicity properties, and reference cases. Implementation
and tests should then be derived from those constraints and reviewed against the mathematics.

The practical rule is: **document intent broadly, specify contracts precisely, generate only
where rules are mechanical, and test every generated or handwritten behavior at its boundary.**
Do not require every note to become a generator input. Start with one small vertical slice and
measure whether generation reduces drift without making the design harder to understand.

Standardized task specification
--------------------------------

Use the repository's maintained
`AI-assisted engineering task specification <https://github.com/eric2003/OneFLOW/blob/hexin/doc/AI_TASK_TEMPLATE.md>`_
for substantial implementation and refactoring work. A
`Chinese version <https://github.com/eric2003/OneFLOW/blob/hexin/doc/AI_TASK_TEMPLATE_zh_CN.md>`_
is available for Chinese-language task preparation.

The specification defines the required project context, observable objective,
whole-system priority assessment, pre-change investigation, implementation
constraints, verification levels, documentation and scientific provenance,
scope control, acceptance criteria, and delivery evidence.

This architecture review remains the system-level decision framework. The task
specification defines how an individual work item is expressed and executed.
Before implementation, read `AGENTS.md`, `doc/AI_CONTEXT.md`, this review,
and the relevant subsystem guide; verify the current `hexin` HEAD. Before
completion, report exact-revision evidence and reassess whether the work changes
the system-wide priority order.

Whole-system architecture map
-----------------------------

The following map is a review checklist for CFD/CAE software. It is intentionally broader
than the currently verified OneFLOW feature set.

.. list-table::
   :header-rows: 1
   :widths: 22 43 35

   * - Concern
     - Questions to answer
     - OneFLOW review entry
   * - Project, case, and configuration
     - What is a project versus a case and a run? Where do parameters come from? How are
       defaults, types, validation, paths, and overrides defined?
     - Runtime guide; R1 case configuration and I/O roots.
   * - Geometry and mesh generation
     - How are geometry entities represented? How are mesh layouts, distributions,
       boundaries, and mesh-generation tasks specified?
     - Mesh architecture and grid-topology guides; mesh layout/workflow code.
   * - Mesh data and topology
     - What is canonical data versus derived connectivity? How are dimensions, cell/face
       relationships, boundary IDs, ghost entities, and cache invalidation defined?
     - Mesh architecture and grid-topology guides; mutable source arrays remain a review
       concern.
   * - Physical models and numerical methods
     - Where are governing equations, constitutive/closure models, boundary conditions,
       discretization, fluxes, limiters, and time integration selected and validated?
     - Trace representative Grid and Solve task paths before proposing broader boundaries.
   * - Linear algebra and solver infrastructure
     - How are matrices/operators, preconditioners, convergence criteria, tolerances, and
       solver configuration owned and tested?
     - Review through concrete solver call paths; do not infer coverage from build flags.
   * - Runtime lifecycle and state
     - Which state is process-scoped, case-scoped, task-scoped, or cached? What happens on
       repeated cases, partial failure, restart, and teardown?
     - SimuContext, SimuImp, global registries, accelerator state, and R1/R3.
   * - Parallel and accelerator execution
     - What is the ownership and synchronization contract across ranks, devices, streams,
       and backend-specific resources? Which test suites actually exercise each backend?
     - Runtime lifecycle plus backend-specific contract and regression evidence.
   * - Input, output, and artifact provenance
     - Where are inputs read and outputs written? Can a consumer identify the exact
       producer artifact? Are source inputs protected from run-time mutation?
     - Runtime I/O roots and Workbench run-location/artifact contracts.
   * - Workflow and orchestration
     - How are independent tasks distinguished from dependent tasks? How are dependencies,
       failure propagation, cancellation, retries, and persisted run results represented?
     - Workbench architecture and R4; planned features must not be described as implemented.
   * - Visualization and data exchange
     - Can runtime data be summarized independently from serialization? Are schemas
       versioned, and can tools inspect data without coupling to the solver?
     - Mesh snapshot, summary, serialization, and visualization adapter work.
   * - Verification and reproducibility
     - Are unit, integration, numerical regression, platform/backend, and performance
       checks distinct? Are baselines, tolerances, and measurement conditions explicit?
     - tests/, test/, CI workflows, baseline documentation, and evidence records.
   * - User and developer experience
     - Can users create/open projects, understand run state, diagnose failures, and inspect
       outputs? Can contributors find build, test, architecture, and recovery guidance?
     - Workbench architecture, published manual, runbooks, and documentation consolidation.
   * - Build, packaging, and extension boundaries
     - Are optional dependencies and backends accurately represented? Can components be
       built and tested independently where useful?
     - Root and component CMake files, ports, and platform-specific CI evidence.

Current cross-system risks to keep visible
------------------------------------------

These are review priorities, not a declaration that every listed area is broken. Recheck
the implementation and exact branch revision before changing status.

1. **Implicit global state across runtime boundaries.** SimuContext provides an explicit
   per-run context, but some task implementations still reach legacy free functions and
   globals. The architectural question is which state should be context-owned and which
   truly belongs to the process; passing a context without using it does not decouple code.
2. **Configuration parsing, validation, and application.** The legacy parameter database
   is shared by multiple subsystems. Parameter syntax, stored type, domain validation, path
   resolution, and mutation of global runtime state are distinct responsibilities. R1 is
   the active investigation; do not generalize a resolver before the current contract is
   covered by tests.
3. **Mesh canonical data and derived topology.** Publicly mutable arrays can make derived
   connectivity and caches stale unless rebuild/invalidation rules are followed. Preserve
   documented topology and boundary invariants while identifying concrete mutation paths
   and their tests.
4. **Workflow artifacts and execution locations.** Source cases, staged run inputs,
   generated meshes, solver results, logs, and workflow-produced artifacts must not be
   conflated. A workflow dependency should consume an identified artifact from its producer,
   not rely on an accidental working directory.
5. **Lifecycle across repeated and failed runs.** Process-level and case-level resources
   have different lifetimes. Changes to setup, accelerator allocation, exception cleanup,
   and multi-case execution need failure-path and repeated-case evidence.
6. **Architecture status versus evidence.** AI handoff pages, roadmaps, and design notes
   reduce rediscovery cost but cannot substitute for focused tests, numerical regression,
   or CI on the exact commit.
7. **Documentation discoverability and authority.** The published manual should teach
   stable concepts and user procedures; engineering notes should retain detailed decisions
   and evidence; AI instructions should point to both without becoming a competing
   architecture source.

Milestone review mechanism
--------------------------

At each meaningful milestone (for example, a completed refactoring slice, a changed
module boundary, a new workflow capability, or a regression investigation), add a short
review to the relevant roadmap or stage record. Do not create a report for every command.

Answer these questions:

#. What changed in behavior, dependency direction, ownership, or testability?
#. Which concern in the whole-system map did it improve?
#. Which concern did it not address, even if adjacent?
#. Did the work reveal a larger issue or change the priority order?
#. What evidence supports the status? Name exact tests and commit-specific CI.
#. What is the next highest-value action, and what would make it complete?

Use consistent status labels:

* **Verified**: supported by source, a named test, or exact-commit CI; cite the evidence.
* **Implemented, unverified**: code exists but required validation is missing.
* **Partial**: some acceptance criteria pass; at least one remains open.
* **Planned**: not implemented.
* **Blocked**: depends on a named missing seam or external condition.
* **Deferred**: intentionally postponed, with the reason.

Architecture progress should be judged by reduced coupling, clearer contracts, reliable
lifecycle behavior, useful tests, reproducible engineering results, and understandable
user workflows—not by the number of classes, interfaces, or documents created.

Relationship to the rest of the handbook
----------------------------------------

* :doc:`runtime_architecture` maps the current execution lifecycle and configuration
  boundaries.
* :doc:`mesh_architecture` and :doc:`grid_topology` describe mesh information and topology
  contracts.
* :doc:`workbench_architecture` describes project, workflow, run, and artifact concepts.
* The AI project context is the concise current-state handoff: https://github.com/eric2003/OneFLOW/blob/hexin/doc/AI_CONTEXT.md
* The runtime refactoring roadmap records staged implementation and acceptance evidence:
  https://github.com/eric2003/OneFLOW/blob/hexin/doc/plans/runtime-refactoring-roadmap.md
* The documentation consolidation plan defines where engineering records and published
  guidance belong:
  https://github.com/eric2003/OneFLOW/blob/hexin/doc/plans/oneflow-documentation-consolidation-plan.md


Code-level quality signals: const contracts and pointer ownership
-----------------------------------------------------------------

Terms
-----

* **Const-correctness** is the practice of expressing read-only access in C++ types.
  A const-qualified member function promises not to modify the observable state of the
  object through that member function; a const reference or pointer prevents mutation
  through that particular access path. This is a compile-time contract, not a guarantee
  that no state anywhere changes.
* **Read-only interface** describes an API surface that permits observation but not
  mutation through that surface. A const view can still refer to an object that another
  alias mutates, so the lifetime and aliasing contract matters.
* **Ownership semantics** describe which object is responsible for keeping a resource
  alive and releasing it. ``std::unique_ptr`` expresses exclusive ownership; a reference
  or raw pointer often expresses non-owning access, but does not encode that intent
  reliably by itself.
* **Borrowed/observer pointer** is a pointer used to access an object owned elsewhere.
  It may be nullable, and its validity depends on the owner's lifetime. A raw pointer is
  not automatically wrong: legacy APIs and graph relationships may need non-owning links.
* **Nullability** asks whether absence is a valid state. If null is meaningful, document
  or model it explicitly. If a value must exist, a reference or a stronger handle may
  express the contract better.
* **Lifetime safety** concerns whether an object remains alive for every use. Smart
  pointers help express ownership, but do not by themselves solve dangling references,
  cycles, concurrent access, or invalidation of views into containers.

Why these are useful signals
----------------------------

Const-correctness can reveal accidental mutation and make dependencies easier to reason
about. Explicit ownership can reveal unclear resource lifetimes, possible leaks, double
deletion, and APIs that blur owning and non-owning roles. These properties are part of
maintainability and testability, but neither can be reduced to a single percentage or
treated as a complete measure of software quality.

A practical inventory should distinguish at least these categories:

.. list-table::
   :header-rows: 1
   :widths: 24 38 38

   * - Signal
     - What to count or inspect
     - Important limitation
   * - Const-correctness
     - Read-only member functions, const input references, const observer pointers/views,
       and mutating APIs with explicit names/contracts.
     - Counting const tokens is misleading; some data is intentionally mutable and
       a const object can still expose mutable aliases.
   * - Ownership clarity
     - Owning unique_ptr/shared_ptr values, documented borrowed pointers, raw
       owning pointers, and APIs that transfer ownership.
     - Smart-pointer counts alone do not reveal cycles, shared-ownership misuse, or
       unclear lifetime dependencies.
   * - Nullability
     - Raw pointers and nullable handles classified by whether null is valid and how
       callers handle it.
     - Replacing every pointer with a reference is wrong when absence is meaningful.
   * - Lifetime and invalidation
     - Views, iterators, references into containers, cached topology, and resources whose
       validity ends after rebuild or teardown.
     - Requires call-path analysis and targeted tests; text search cannot prove safety.
   * - Evidence
     - Focused contract tests, sanitizer/static-analysis findings, compiler diagnostics,
       and exact-revision CI.
     - Passing tests demonstrate covered behavior, not all possible lifetime correctness.

Suggested reporting
-------------------

Track both a **coverage inventory** and a **risk inventory**. For example, report the
number of audited APIs, how many have explicit ownership/nullability contracts, how many
read-only operations are const-qualified where appropriate, and how many high-risk raw
owning pointers remain. Include the audited directories, exclusions, and review date.
Treat unknown cases as unknown rather than silently counting them as safe.

A percentage is meaningful only when its denominator is defined. Prefer metrics such as
"borrowed pointers reviewed in codes/geometry" or "public APIs with documented ownership
and nullability" over "percentage of the whole codebase that is safe". Track high-risk
unresolved cases separately, because one raw owning pointer on a critical failure path may
matter more than many low-risk declarations.

For OneFLOW, apply this incrementally to real call paths. Start with geometry code that
constructs and consumes interface mappings, then expand to mesh topology and runtime
lifecycle when the first inventory shows a useful pattern. Do not launch a repository-wide
mechanical pointer conversion. Each change should preserve numerical behavior, state
which alias is allowed to mutate the object, identify the owner, and add a test when the
contract can be observed.


A broader quality model for refactoring
----------------------------------------

Const-correctness and ownership clarity are useful code-level signals, but they cover
only a small part of software quality. ISO/IEC 25010:2023 defines a broader product
quality model; its characteristics are a useful vocabulary, not a demand that every
project metric be collected. For OneFLOW, the model should be adapted to scientific
computing and to the current subsystem being changed.

.. list-table::
   :header-rows: 1
   :widths: 20 43 37

   * - Quality dimension
     - Practical evidence or metric
     - How to avoid misleading results
   * - Functional correctness
     - Contract tests, known-answer tests, integration tests, defect recurrence tests.
     - A test count or pass rate says nothing about untested behavior.
   * - Numerical validity
     - Reference-case errors, conservation residuals, convergence behavior, mesh/topology
       invariants, and changes against an accepted baseline.
     - Record mesh, parameters, compiler/backend, tolerances, and comparison method.
   * - Reliability and lifecycle
     - Repeated-run tests, setup/teardown tests, failure-path tests, resource diagnostics,
       and recovery behavior.
     - Include normal reuse and realistic failures; do not optimize for fantastical cases.
   * - Modularity and coupling
     - Dependency direction, cross-layer includes/calls, global-state reachability, and
       whether a component can be tested without starting the whole solver.
     - Fewer dependencies are not automatically better if the design becomes indirect.
   * - Cohesion and responsibility
     - Whether a module has one understandable reason to change, and whether algorithms,
       data ownership, I/O, and orchestration are mixed together.
     - Avoid treating class size or file length as a verdict by itself.
   * - API and state contracts
     - Const-correctness, ownership, nullability, invalidation rules, units, valid states,
       error reporting, and preconditions/postconditions.
     - Review contracts at call sites; declarations alone may not reveal aliasing or
       lifecycle assumptions.
   * - Testability and verification
     - Important contracts with focused tests, registered test discovery, integration
       coverage, numerical regression coverage, and actionable failure messages.
     - Coverage percentage is a map of exercised code, not proof of correctness.
   * - Complexity and change risk
     - Cyclomatic/cognitive complexity for critical routines, duplicated policy, branching,
       and the number of unrelated components touched by a change.
     - Use trends and review hotspots; avoid arbitrary repository-wide thresholds.
   * - Performance and resource use
     - Profiled runtime, memory high-water mark, allocation count, scaling, and hotspot
       benchmarks on representative cases.
     - Compare equivalent workloads and configurations; benchmark only where performance
       matters to users.
   * - Portability and buildability
     - Supported compiler/platform builds, clean configure/build, registered tests, and
       dependency configuration.
     - A successful build on one machine does not establish cross-platform compatibility.
   * - Reproducibility and provenance
     - Recorded inputs, executable/build identity, case location, output artifacts, and
       stable regression comparisons.
     - Distinguish deterministic requirements from expected floating-point variation.
   * - Usability and diagnosability
     - Clear configuration errors, useful logs, discoverable workflows, actionable failure
       messages, and the time needed to reproduce a reported issue.
     - Evaluate actual user tasks, not the number of UI controls or documentation pages.
   * - Documentation quality
     - Documented public contracts, algorithm rationale, examples that are tested where
       practical, and links from the developer entry point to authoritative subsystem docs.
     - Documentation volume is not a quality metric; correctness and discoverability matter.

Metric selection and use
------------------------

Do not turn this table into a large scorecard that every change must fill out. For each
refactoring slice, select a small set of signals based on its risk:

* **Every slice:** behavior preserved or intentionally changed; focused tests; relevant
  integration or regression evidence; exact revision and test registration.
* **A boundary or ownership change:** dependency direction, owner/borrower relationship,
  lifetime and invalidation contract, and whether the component can be exercised in
  isolation.
* **A numerical algorithm change:** reference outputs, tolerances, convergence or
  conservation evidence, and a representative performance check if the algorithm is
  performance-sensitive.
* **A workflow or runtime change:** repeatability, failure cleanup, artifact provenance,
  logs, and a realistic end-to-end path.
* **A performance change:** benchmark method, workload, environment, and before/after
  result; do not infer speed from cleaner-looking code.

For each chosen metric, define its denominator and collection method before comparing
numbers. Keep **leading indicators** (for example, unclear ownership contracts or
unregistered tests) separate from **outcomes** (for example, numerical regressions,
crashes, or measured runtime). A rise in test count, const usage, or documentation count
is not itself proof of better software.

For the current geometry work, the useful combination is not a single quality score. It
is: preserve the point-matching algorithm's observable behavior; state its tolerance and
candidate-search contracts; ensure tests are discoverable; understand how its borrowed
grid views and owned ADT interact; and only then decide whether a larger dependency or
lifecycle boundary needs redesign. This connects a small local change to testability,
numerical correctness, and maintainability without expanding it into a repository-wide
cleanup.


Point locator: index and coordinate-store lifecycle
--------------------------------------------------

The point locator maintains two representations of inserted points: the ADT spatial
index stores point nodes and IDs, while ``xCoor``, ``yCoor``, and ``zCoor`` store the
coordinates used by point accessors. These structures are one logical data set even
though they are represented separately.

The query algorithm has two stages:

#. The ADT broad phase retrieves candidates inside an axis-aligned box centered on the
   query point and expanded by the active tolerance.
#. The narrow phase computes squared Euclidean distance and accepts only candidates with
   ``distance_squared < tolerance * tolerance``. It returns the closest accepted point;
   equal-distance ties depend on candidate traversal order.

The squared-distance comparison avoids a square root for every candidate. The strict
inequality is observable behavior at the tolerance boundary and should not be changed
without an explicit compatibility decision.

A lifecycle invariant follows from the two representations: whenever initialization
replaces the ADT, the coordinate arrays for points inserted through this locator must
also be cleared. Otherwise the new index can be empty while the old coordinate arrays
remain populated; the next insertion can receive an ID based on stale coordinates, and
the index and coordinate store no longer describe the same set of points.

The reinitialization test therefore checks more than safe destruction: after rebuilding,
the point count is zero, a previously inserted point is no longer found, and the next
insertion starts again at ID zero. This is a targeted example of a useful quality
measure: **cross-representation invariant coverage**. It is more meaningful here than
counting lifecycle tests or smart pointers.

This change intentionally does not alter tolerance policy, tie-breaking behavior,
grid bounding-box calculation, or the ADT search algorithm. Those are separate
contracts and should be evaluated independently. The locator's inserted-point store is
also distinct from the source grid nodes used to establish initial search bounds; clearing
the former does not clear or modify the source grid.


Special tolerance and broad-phase bounds
----------------------------------------

``InitializeSpecial(grid, toleranceIn)`` must keep the broad and narrow phases in
agreement. The narrow phase accepts a point only when its squared distance is strictly
less than ``toleranceIn * toleranceIn``. The broad phase first asks the ADT for points in
a box expanded by that same tolerance. If the ADT root bounds are expanded only by the
grid-derived default tolerance (``minimum_spacing / 4``) while the query uses a larger
explicit tolerance, a point can satisfy the exact distance test but never reach the
narrow phase because it was excluded by the index bounds.

The implementation now calculates the source grid bounding box and rebuilds the ADT
through the explicit-bounds initialization path with ``toleranceIn``. This keeps the
index bounds, query window, and exact distance threshold aligned, and also reuses the
existing initialization path that resets the inserted-point coordinate store. The normal
single-grid and multi-grid tolerance policies remain unchanged.

This is an algorithm-contract fix, not a general redesign of the locator. It illustrates
why a spatial index has two correctness conditions: its candidate structure must not omit
any point that can pass the exact predicate, and the exact predicate must preserve the
intended tolerance semantics. The registered ``SpecialToleranceExpandsAdtBounds`` lifecycle test exercises this
contract with a synthetic grid: it inserts a point outside the raw grid bounds but inside
the explicitly expanded ADT bounds, then verifies that a nearby query finds it.
