OneFLOW Runtime Architecture and Development Guide
===================================================

Purpose and status
------------------

This page is the system-level map for developers who need to understand the intent,
architecture, configuration boundaries, program entry points, and execution flow before
changing individual classes. It is based on the `hexin` source tree and the curated mesh
and Workbench design pages.

This is an initial architecture map, not a claim that every intended boundary is fully
enforced. Where current code still relies on global state, mutable arrays, or legacy
configuration syntax, that limitation is called out explicitly. Recheck the named code
paths when implementation changes.

Design intent
-------------

OneFLOW combines a C++ simulation runtime with mesh generation and conversion, case-based
configuration, numerical execution, regression tooling, and an optional Python/PySide6
Workbench. The intended boundary is:

* The C++ runtime owns simulation and numerical execution.
* A case directory supplies the control inputs and other files required by that run.
* Mesh and solver modules implement domain behavior; command-line and Workbench tools
  select cases, prepare execution locations, and inspect artifacts.
* Tests state behavior contracts at the smallest useful boundary, while end-to-end
  regressions check that the assembled application still behaves correctly.
* Documentation records design intent and current contracts separately from proposed
  future capabilities.

The Workbench is an orchestration and inspection layer, not a second solver. Mesh summary
and serialization code should describe runtime data without making the core mesh model
depend on JSON or GUI libraries.

System map
----------

.. code-block:: text

   User / CLI / Workbench
            |
            v
   Case selection and run preparation
            |
            v
   OneFLOW executable: main()
            |
            v
   Simulation::Run()
       |                         |
       | no case-directory args  | case-directory args
       v                         v
   lightweight test selector   parse command-line options
                                 |
                                 v
                           process environment
                                 |
                                 v
                         one SimuImp per case
                                 |
                                 v
                         SimuContext / case setup
                                 |
                                 v
                         control task resolution
                                 |
                                 v
                         TaskRegistry -> ISimuTask
                                 |
                                 v
                         mesh / solver / post task
                                 |
                                 v
                         case teardown and process finalization

This is a conceptual flow, not a complete call graph. Specific task bodies may still call
legacy free functions and access shared runtime state. See the code entry points below.

Code entry points and responsibilities
--------------------------------------

Application entry:

* `codes/main/src/OneFlow.cpp` contains the small `main(argc, argv)` entry point. It
  constructs `Simulation` and returns `Simulation::Run()`.
* `codes/main/src/Simulation.cpp` processes command-line arguments and chooses between
  the lightweight test-selector path and full case execution.
* `codes/main/src/SimuImp.cpp` owns a case lifecycle: setup, task resolution, execution,
  and teardown. Its exception path attempts teardown before rethrowing the original error.
* `codes/main/src/SimuContext.cpp` contains testable context state and selected state
  operations. Environment setup is separated into `SimuContextEnv.cpp`.
* `codes/main/src/SimuTaskReg.cpp` registers task implementations. Concrete task classes
  adapt existing operations to the `ISimuTask` interface.
* `codes/io/src/ConfigLoader.cpp` parses the legacy parameter syntax into a
  `ConfigDocument`. Variable-sized arrays require an injected size resolver when their
  size is not a numeric literal.

The build is configured from the root `CMakeLists.txt`, which selects the third-party,
runtime, test, and GUI subdirectories. `codes/CMakeLists.txt` collects runtime sources
and defines build options for optional acceleration and parallel facilities. Build-time
options do not prove that a backend is available or validated at runtime; use the build
summary and relevant tests as evidence.

Command-line execution modes
----------------------------

The current `Simulation::RunImpl()` has three argument-count paths:

* With only the executable argument, it enters the lightweight test selector and chooses
  a registered test, with `ONEFLOW_DEFAULT_TEST` as an optional environment override.
* With exactly one additional argument, it reports an argument-count error.
* With two or more additional arguments, it parses the case-directory arguments, performs
  process-level initialization, and runs each resolved case directory through its own
  `SimuImp` instance.

For example, `OneFLOW.exe r <caseDir>` has three process arguments and follows the full
simulation path. Case directories passed to one process are executed in command-line order; every argument entry is executed, so repeating the same path runs that case again. Each entry gets a separate ``SimuImp`` / ``SimuContext`` and case teardown, while process runtime (including MPI initialization) is initialized once and finalized after the sequence. ``SimuImp::RunCase()`` can be called repeatedly inside that process lifetime; it owns case teardown, not MPI lifecycle. ``SimuImp::Run()`` is the single-case compatibility entry point that owns one process lifecycle. If process initialization or case setup/run throws, the owning boundary attempts rollback/finalization before propagating the original failure. The compatibility ``SimuContext::SetupEnvironment()`` path rolls back a failed case setup, and ``TeardownEnvironment()`` finalizes process runtime even if case teardown throws. A failed case stops the sequence; later entries are not run. This is not parallel multi-case execution and does not itself express artifact dependencies. Thirteen focused lifecycle tests cover ordering, repeated paths, rollback only on setup failure, exception preservation, and stop-on-failure behavior. ``SimuContextEnv.cpp`` compiled separately with MSVC. The full suite and executable link remain unverified. Use the Workbench architecture page for the distinction between sequential multi-case runs and dependency workflows.

The process environment and each case environment have different lifetimes. The
`Simulation` path initializes process-level state once, then constructs a case runner per
directory. `SimuImp::RunCase()` owns case teardown, including exceptional exits. Changes
to these lifetimes should be protected by tests that cover both successful execution and
failure during setup or task execution.

Configuration boundaries
-------------------------

OneFLOW configuration is not a single undifferentiated settings file. Keep these layers
separate when changing configuration behavior.

Build configuration
   CMake cache options, compiler selection, dependency locations, and optional backend
   switches determine what is compiled and linked. These settings belong to configure and
   build operations, not to a particular case's physical model.

Case runtime configuration
   Case control files and legacy parameter files select runtime tasks and numerical,
   geometry, or mesh-generation parameters. The task name is resolved from the control
   database by `SimuContext::ResolveTaskFromControl()`; current task names include
   `Solve`, `Grid`, `WallDist`, `Partition`, `FunctionTest`, `Theory`,
   `ToyModel`, and `PostTask`. `ConfigLoader` parses legacy text syntax into a
   `ConfigDocument`; it is not itself the owner of all runtime meaning. The code that
   resolves a parameter must define its type, default, valid range, and effect.

   Case I/O roots are applied during `SimuContext::SetupCaseEnvironment()`. The current
   configuration keys include `projectRoot`, `inputRoot`, `outputRoot`,
   `meshInputRoot`, `flowInputRoot`, `meshOutputRoot`, `flowOutputRoot`, and
   `ioPathInfoFileName`. If `meshInputRoot` is empty or `default`, a non-empty
   `sourceGridCaseDir` is used as the legacy mesh-input fallback. Generic input/output
   roots are applied before the mesh/flow-specific roots. Missing optional root settings
   retain the existing case-relative defaults.

   Optional path strings are read through
   `SimuContextDetail::ReadOptionalStringParameter()` in
   `codes/main/include/SimuContextParameter.h`. The helper returns an empty string for a
   missing entry, verifies that a present entry is `HX_STRING`, and throws
   `std::invalid_argument` naming the parameter when its declared type is wrong. This
   narrows the earlier unchecked typed-read behavior for these optional path settings;
   it does not make every control parameter type-safe, nor does the helper test validate
   global database binding or the complete case-setup lifecycle.

   The legacy mesh input fallback is split into two explicit decisions.
   `ShouldUseLegacySourceGridCaseDir()` decides whether the configured
   `meshInputRoot` permits fallback lookup (only empty or `default`). The
   `sourceGridCaseDir` parameter is read lazily only in that case, then
   `ApplyLegacySourceGridCaseDir()` chooses the effective value: a non-empty,
   non-`default` legacy value replaces the fallback token; otherwise the original
   mesh input value is preserved. Explicit roots, including `case` and ordinary
   relative paths, therefore do not consult or get overridden by the legacy setting.
   Focused tests cover both the lookup decision and the effective-value contract.

   Existing `Prj` tests cover generic/domain root resolution, inheritance, diagnostics,
   and resetting roots when binding another case. The helper-level
   `simu_context_parameter_test` covers missing, valid-string, and invalid-type lookup.
   A direct `SetupCaseEnvironment()` integration test is still needed to prove the
   precedence between explicit `meshInputRoot` and `sourceGridCaseDir`, as well as the
   interaction with control-file loading and sequential case setup. Keep the documented
   contract aligned with `SimuContextEnv.cpp` and the path-resolution tests.

Workbench project and workflow configuration
   Project definitions describe an engineering workspace; workflow definitions describe
   execution steps and their dependencies; run records describe individual attempts and
   their results. Definitions and run history must remain separate. A case catalog helps
   discover cases but must not silently determine the actual run directory.

Path and artifact semantics
   Source case inputs, staged run inputs, generated mesh files, solver results, and logs
   have different roles. Any operation that can modify files should make its execution
   location explicit and should isolate source inputs when that operation's contract
   requires a staged copy. Producer/consumer workflows should resolve declared artifacts,
   not depend on an accidental current working directory.

Patterns already present
-----------------------

Use pattern names only where they describe real code, not as a goal to apply everywhere.

Command / Strategy boundary
   `ISimuTask` represents an executable task, and concrete tasks adapt existing operations
   to that interface. This gives task selection a stable seam, but some implementations
   still delegate to legacy free functions and global state.

Registry / Factory
   `TaskRegistry` maps control-file task names to creators. The registry provides a
   centralized construction point. Registration behavior, duplicate-name policy, and
   initialization order should remain explicit and testable.

Context object
   `SimuContext` groups selected run state and provides an injectable seam for context
   data. Its value depends on progressively moving relevant state out of implicit globals;
   adding a context parameter without using it does not by itself decouple a task.

Adapter
   `ConfigLoader` translates legacy text syntax into `ConfigDocument`. Similar adapter
   boundaries are appropriate for runtime mesh snapshots and serialization formats.
   Parsing syntax and interpreting domain semantics should remain separate responsibilities.

Explicit lifecycle / RAII
   `Simulation` and `SimuImp` use `std::unique_ptr` for owned objects, and case teardown
   is expressed as a lifecycle operation. Continue to make ownership and cleanup visible;
   avoid adding wrapper classes that do not clarify a real lifetime or responsibility.

These are existing patterns with uneven adoption, not a claim that the whole codebase is
already dependency-inverted or free of global state.

Architecture rules for ongoing refactoring
-------------------------------------------

1. Follow the real call path before changing an interface. Inspect direct callers, ownership,
   initialization order, and teardown behavior.
2. Establish the behavior contract before extracting a new class or renaming a module.
   Tests should distinguish intended behavior from plausible but incorrect behavior.
3. Prefer value snapshots or narrow read-only interfaces for data that crosses a module
   boundary. Do not expose mutable implementation containers merely to simplify a test.
4. Keep canonical state distinct from derived data. Any derived adjacency, cache, or mapping
   needs an explicit rebuild or invalidation contract.
5. Keep numerical algorithms in the runtime; orchestration code should select and launch
   work rather than reimplement solver logic.
6. Separate parsing from validation and semantic interpretation. Report invalid or unsupported
   input explicitly instead of silently assigning an unrelated default.
7. Refactor one independently verifiable boundary at a time. Avoid broad rewrites without a
   concrete behavior or maintainability problem.
8. Keep source code ASCII-only, preserve compatibility where required, and validate the
   relevant Windows and Linux paths when the change affects them.

Test design as an architecture probe
------------------------------------

A test that requires the whole solver to validate a small rule may indicate an overly broad
dependency boundary, hidden global state, or setup responsibilities that are hard to isolate.
It is a signal to investigate, not automatic proof that a new abstraction is needed.

Use a layered test ladder:

* Unit tests for deterministic parsing, mapping, lifecycle, and invariant rules.
* Component tests for collaboration between context, task registry, mesh services, and
  configuration adapters.
* Integration tests for real case setup, generated artifacts, and runtime paths.
* CFD/grid regression tests for end-to-end behavior and numerical results.

Tests should assert externally meaningful contracts rather than mirror private data
structures. If a test is difficult to write, first identify the source of the difficulty:
algorithmic complexity, hidden dependencies, broad setup, or a genuinely necessary system
integration.

Documentation change loop
-------------------------

For each meaningful refactor:

#. State the observed problem and the intended behavior.
#. Add or update the smallest test that distinguishes correct behavior from the defect.
#. Inspect callers and ownership, then make the smallest coherent code change.
#. Update the relevant architecture page when responsibilities, lifecycle, data contracts,
   configuration, entry points, or user-visible behavior change.
#. Record an Architecture Decision Record only when the choice has durable consequences or
   meaningful alternatives. Do not create an ADR for every local fix.
#. Run the relevant tests and record exactly what was executed. Distinguish local results,
   Linux CI, Windows CI, and numerical regression evidence.
#. Update user instructions only when supported behavior or operating procedures change.

The code, tests, and documentation should be reviewed together. A document may describe an
intended boundary before it is fully implemented only if it labels that boundary as
proposed and names the gap.

Related design pages
--------------------

* :doc:`mesh_architecture` — mesh snapshots, canonical data, and serialization boundaries.
* :doc:`grid_topology` — topology invariants, boundary normalization, and derived connectivity.
* :doc:`workbench_architecture` — projects, workflows, run locations, artifacts, and process lifecycle.
* :doc:`info` — project scope and documentation entry points.
* `Documentation consolidation plan <https://github.com/eric2003/OneFLOW/blob/hexin/doc/plans/oneflow-documentation-consolidation-plan.md>`_ — document ownership and migration status.
* `Runtime refactoring roadmap <https://github.com/eric2003/OneFLOW/blob/hexin/doc/plans/runtime-refactoring-roadmap.md>`_ — staged refactoring goals and acceptance criteria.
