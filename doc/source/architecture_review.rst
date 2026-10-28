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

A reusable prompt for AI-assisted engineering
---------------------------------------------

Use this prompt at the start of a substantial task and again before declaring it done:

.. code-block:: text

   Work from the latest hexin HEAD unless the task explicitly names another branch.
   First read AGENTS.md, doc/AI_CONTEXT.md, this architecture review, and the relevant
   subsystem guide. Confirm the current branch and revision.

   Before editing:
   - State the user-visible or maintenance problem in concrete terms.
   - Map the actual entry point, direct callers, state/data ownership, lifecycle, failure
     cleanup, tests, and downstream consumers.
   - Locate the issue in the whole-system architecture map below.
   - Separate verified behavior, suspected gaps, and proposed future design.
   - Identify the smallest testable change and the evidence needed to accept it.
   - Check whether the issue is local or a symptom of a larger unresolved boundary.

   While working:
   - Preserve compatibility unless a behavior change is explicit and tested.
   - Do not introduce abstractions without a concrete consumer and benefit.
   - Keep C++ source ASCII-only.
   - Do not claim builds, tests, or CI results that were not observed for the exact revision.

   Before finishing:
   - Report what changed and what deliberately did not change.
   - Report test commands, results, failures, and exact-revision CI evidence.
   - Revisit the subsystem and whole-system maps: what architectural risk was reduced,
     what remains open, and should priorities change?
   - Update the relevant architecture page and roadmap, and update doc/AI_CONTEXT.md
     only when the overall status or next action materially changes.

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
