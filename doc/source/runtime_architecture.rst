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

System-level design principles
-------------------------------

The runtime should be understood through a small set of boundaries before examining a
particular class or function:

1. **Numerical ownership:** the C++ runtime owns numerical state and computation. The
   CLI and Workbench select, stage, launch, and inspect runs; they do not reimplement
   solver algorithms.
2. **Lifecycle ownership:** process setup/finalization and per-case setup/teardown have
   different lifetimes. A case may fail without silently skipping cleanup of the process
   resources it owns.
3. **Data ownership:** live mesh and solver objects are operational state. A snapshot is
   a stable read model; JSON, VTU, or other output formats are adapters, not the core
   data model.
4. **Dependency direction:** domain modules should not need a GUI or a particular output
   serialization library to perform numerical work. Presentation and orchestration depend
   on domain contracts, not the reverse.
5. **Evidence at boundaries:** unit tests pin small contracts; integration tests check
   module interaction; grid and CFD regressions check behavior in representative cases.
   A passing test count is useful evidence only when its coverage and assertion contract
   are understood.

These principles form the review checklist for local refactoring. A change should state
which boundary it clarifies, which behavior must remain unchanged, and what evidence
would detect a regression. Avoid creating a new abstraction when a simpler name, explicit
lifecycle, or focused contract test already makes the responsibility clear.

Source comments and generated documentation
-----------------------------------------

Use source comments for contracts needed at declarations and call sites: purpose,
units, ownership, valid ranges, side effects, failure behavior, and links to the
design page. Use Doxygen-compatible comments for stable public interfaces and
algorithm entry points; avoid duplicating long tutorials in headers.

Keep ``doc/source`` as the canonical home for system architecture, design rationale,
equations and symbol definitions, assumptions, diagrams, references, and validation
evidence. Generated API pages complement these pages; they do not replace them.

A staged route is preferable to immediately wiring a large documentation toolchain:

1. Standardize Doxygen-compatible comments on one stable module or interface.
2. Generate API documentation into a build directory and check that generation is
   deterministic and does not dirty the source tree.
3. If API pages need to appear inside Sphinx, evaluate Doxygen XML plus Breathe only
   after the initial output is useful and repeatable.
4. Keep formulas in Sphinx math directives, and use small source-controlled flow
   diagrams to show lifecycle, data flow, and ownership. Each diagram should distinguish
   current behavior from proposed design and link to implementation entry points.
5. Treat generated documentation in CI as a reproducible artifact; do not commit
   generated output unless the project explicitly adopts that publication model.

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

The simulation-specific schema and invocation policy are implemented in ``codes/io/src/SimulationCommandLine.cpp`` and declared in ``codes/io/include/SimulationCommandLine.h``. The generic parser has no dependency on ``Prj`` or solver runtime types. ``PrjCmdLineArgs.cpp`` is limited to converting canonical options to the legacy ``CmdLineOptions`` shape. The parser core is ``ONEFLOW::CliParser`` in ``codes/basic/include/CliParser.h``. It uses a declared positional schema with minimum/maximum counts, returns named parsed values, preserves repeated values, and generates usage text for validation errors. ``ParseSimulationCommandLineOptions()`` converts the parsed values into a ``SimulationCommandLineOptions`` value object with a single canonical case-directory list. ``Simulation`` and ``SimuContext`` now use this API. ``Prj::ParseCmdLineArgs()`` remains a compatibility adapter that fills the older ``CmdLineOptions`` shape. Named switches, subcommands, and type conversion remain out of scope until an actual command needs them.

.. code-block:: cpp

   CliParser parser( "OneFLOW" );
   parser.AddPositional( "mode", "Run mode" );
   parser.AddPositional(
       "caseDir", "Case directory", 1, CliParser::Unlimited );
   const ParsedCommandLine parsed = parser.Parse( { "r", "caseA", "caseB" } );

``Parse()`` receives arguments after the executable name. The parser preserves the existing mode rule (``d`` enables debug; other values retain run mode). The compatibility adapter still fills the old ``Prj::ParseCmdLineArgs()`` result shape. The small parser intentionally does not yet provide named options or subcommands.

`Simulation::RunImpl()` selects its top-level path through the pure `ClassifySimulationInvocation()` function. One argv entry selects the lightweight test path, two entries preserve the existing argument-count error, and three or more select full simulation. The classifier does not inspect argument text; parsing and validation of the full-simulation arguments remain in `ParseSimulationCommandLineOptions()`. Tests in `tests/main/prj_cmdline_test.cpp` pin these three behaviors so future CLI syntax work can change them deliberately.

The current `Simulation::RunImpl()` has three argument-count paths:

* With only the executable argument, it enters the lightweight test selector and chooses
  a registered test, with `ONEFLOW_DEFAULT_TEST` as an optional environment override.
* With exactly one additional argument, it reports an argument-count error.
* With two or more additional arguments, it parses the case-directory arguments, performs
  process-level initialization, and runs each resolved case directory through its own
  `SimuImp` instance.

For example, `OneFLOW.exe r <caseDir>` has three process arguments and follows the full
simulation path. Case directories passed to one process are executed in command-line order; every argument entry is executed, so repeating the same path runs that case again. Each entry gets a separate ``SimuImp`` / ``SimuContext`` and case teardown, while process runtime (including MPI initialization) is initialized once and finalized after the sequence. ``SimuImp::RunCase()`` can be called repeatedly inside that process lifetime; it owns case teardown, not MPI lifecycle. ``SimuImp::Run()`` is the single-case compatibility entry point that owns one process lifecycle. If process initialization or case setup/run throws, the owning boundary attempts rollback/finalization before propagating the original failure. The compatibility ``SimuContext::SetupEnvironment()`` path rolls back a failed case setup, and ``TeardownEnvironment()`` finalizes process runtime even if case teardown throws. A failed case stops the sequence; later entries are not run. This is not parallel multi-case execution and does not itself express artifact dependencies. Thirteen focused lifecycle tests cover ordering, repeated paths, rollback only on setup failure, exception preservation, and stop-on-failure behavior. ``SimuContextEnv.cpp`` compiled separately with MSVC. The full suite and executable link remain unverified. Use the Workbench architecture page for the distinction between sequential multi-case runs and dependency workflows.

Context ownership transfer
^^^^^^^^^^^^^^^^^^^^^^^^^^^

``SimuContext`` is movable but not move-assignable: moving transfers the case
arguments, task metadata, accelerator-state registry, readiness state, and process
lifecycle marker to the destination. The source is reset to an uninitialized state.
This operation is declared ``noexcept`` because the process lifecycle marker is an
ownership token; once it moves, a later throwing member move could otherwise leave the
source unable to finalize the process runtime while destination construction fails.
The compile-time move-trait assertion in ``tests/main/simu_context_test.cpp`` makes
that requirement executable. If a future context member cannot be moved without
throwing, its ownership design must be reconsidered rather than silently weakening
this guarantee.

The process environment and each case environment have different lifetimes. The
`Simulation` path initializes process-level state once, then constructs a case runner per
directory. `SimuImp::RunCase()` owns case teardown, including exceptional exits. Changes
to these lifetimes should be protected by tests that cover both successful execution and
failure during setup or task execution.

Case lifecycle: success and failure paths
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The current implementation has two nested cleanup responsibilities. The process owner
initializes and finalizes shared runtime resources; the case owner tears down case-scoped
solver, field, parameter, log, and project-path state. The lifecycle helpers in
``codes/main/include/SimuCaseLifecycle.h`` preserve the original execution exception
when cleanup also fails.

.. code-block:: text

   Simulation / SimuImp::Run
       |
       +-- SetupProcessEnvironment
       |      parallel setup -> accelerator runtime setup
       |
       +-- RunCase
       |      |
       |      +-- PreProcess -> InitSimu -> SetupCaseEnvironment
       |      |      case path/log binding
       |      |      -> ReadControlInfo(caseDir)
       |      |      -> read and validate CaseIOSettings
       |      |      -> apply Prj roots -> write io_paths.txt
       |      |
       |      +-- MainProcess -> resolve task -> execute task
       |      |
       |      +-- PostProcess -> TeardownCase
       |
       +-- FinalizeEnvironment
              accelerator finalization -> parallel finalization

If setup or task execution throws inside ``RunCase()``, its cleanup callback attempts
``PostProcess()`` before the original exception is rethrown. The outer process lifecycle
then attempts process finalization. Cleanup failures do not replace an already active
execution failure. If execution succeeds, a cleanup failure is allowed to propagate.
These are observable contracts and should remain explicit when the orchestration boundary
is extracted.

This diagram describes the existing path; it does not imply that case setup is
transactional. In particular, project and log bindings are established before the real
control reader runs. A failed read can therefore leave partial state until the owning
cleanup boundary runs. Direct callers of ``SetupCaseEnvironment()`` do not get the
same automatic rollback that callers of ``SetupEnvironment()`` or ``SimuImp::RunCase()``
receive. Integration tests must invoke the owning production boundary for cleanup claims.

The next R1 acceptance evidence should cover distinct layers:

* A minimal temporary multi-file case exercises the real ``ReadControlInfo(caseDir)``
  and verifies typed database values plus ``log/database.log``.
* A setup-only smoke test uses a verified case under ``test/`` and checks the resulting
  effective roots and path diagnostics without running the solver task.
* Two cases run sequentially in one process to show that the second case does not inherit
  roots from the first.
* A missing or malformed control input verifies cleanup and preservation of the original
  exception.

The existing lifecycle unit tests cover the generic callback ordering and exception
policy, but they do not prove these production integration properties. R1 remains partial
until the production reader, effective path binding, sequential isolation, and failure
cleanup have direct evidence.

Configuration boundaries
-------------------------

Case I/O configuration snapshot
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The case setup boundary reads optional path settings into one ``CaseIOSettings``
value before applying configured roots to legacy ``Prj`` state. The snapshot
includes ``ioPathInfoFileName``, which selects the path-diagnostic output name.
A present setting must be a string; an absent setting uses the historical
empty-value default. Legacy ``sourceGridCaseDir`` is consulted only when
``meshInputRoot`` is empty or ``default``.

The intended ordering is:

.. code-block:: text

   ReadControlInfo(caseDir)
             |
             v
   ReadCaseIOSettings(DataPara)
   [read + type-check all applicable settings]
             |
             v
   ApplyCaseIOSettings(CaseIOSettings)
             |
             +--> Apply project and generic roots
             |
             +--> Apply mesh/flow-specific roots
             |
             v
   Resolve and write io_paths.txt

The declaration for ``SimuContextDetail::ApplyCaseIOSettings()`` is in
``codes/main/include/SimuContextIO.h``; its implementation is compiled from
``codes/main/src/SimuContextIO.cpp``. The implementation applies the validated
snapshot through the existing ``Prj`` API and writes diagnostics to the mesh output root and,
when distinct, the flow output root. The focused test target links the real
``Prj`` path implementation and checks file placement and the effective roots
serialized in each diagnostic file.

This separates parameter interpretation from mutation of global path state. It
does not make setup fully transactional: path application and filesystem output
can still fail, so the surrounding case lifecycle must retain its cleanup path.
These tests pin the post-reader operation; they do not yet prove production
control-reader integration, sequential-case isolation, or failure cleanup.

Root resolution has one important compatibility distinction. ``projectRoot`` is
configured through ``Prj::SetProjectBaseDir()``: an absolute value is normalized
directly, while a relative value is resolved against ``Prj::current_dir`` (the
process working directory captured by project initialization). In contrast,
ordinary relative ``inputRoot``, ``outputRoot``, and domain-specific roots resolve
against the active case directory. The tokens ``case`` and ``default`` select
the case root; ``project`` selects the configured project root and fails if no
project root has been configured. This distinction is implemented in
``codes/io/src/Prj.cpp`` and should not be changed implicitly as part of setup
refactoring. A focused ``PrjIOPaths`` contract test pins the relative
``projectRoot`` behavior separately from case-relative generic output.


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
   `std::invalid_argument` naming the parameter when its declared type is wrong.

   `ReadCaseIOSettings()` snapshots the generic and mesh/flow-specific root settings
   into a value object before the runtime applies them through `Prj`. This separates typed
   parameter lookup and legacy fallback selection from mutation of the process-global path
   fields. The parameter adapter does not normalize paths or apply roots. The separate
   `SimuContextDetail::ApplyCaseIOSettings()` is declared in
   `codes/main/include/SimuContextIO.h` and implemented in
   `codes/main/src/SimuContextIO.cpp`. Keeping path mutation in a source file prevents
   every includer from inheriting the implementation's `Prj` and filesystem dependencies.
   It applies the snapshot and writes path diagnostics; it does not replace the production
   `ReadControlInfo(caseDir)` call.

   The legacy mesh input fallback is split into two explicit decisions.
   `ShouldUseLegacySourceGridCaseDir()` decides whether the configured
   `meshInputRoot` permits fallback lookup (only empty or `default`). The
   `sourceGridCaseDir` parameter is read lazily only in that case, then
   `ApplyLegacySourceGridCaseDir()` chooses the effective value: a non-empty,
   non-`default` legacy value replaces the fallback token; otherwise the original
   mesh input value is preserved. Explicit roots, including `case` and ordinary
   relative paths, therefore do not consult or get overridden by the legacy setting.
   Focused tests cover the full settings snapshot, legacy fallback selection, and the
   fact that an explicit mesh root does not read an invalid legacy parameter.

   Existing `Prj` tests cover generic/domain root resolution, inheritance, diagnostics,
   and resetting roots when binding another case. The helper-level
   `simu_context_parameter_test` covers the complete settings object, missing/valid/invalid
   parameter types, and application of mesh/flow roots to real diagnostic files. A direct
   `SetupCaseEnvironment()` integration test is still needed to verify the real control
   reader, sequential case setup, and failure cleanup together.
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
------------------------

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

Numerical theory, references, and verification evidence
-------------------------------------------------------

Numerical documentation should connect the mathematical model to the implementation
and to observable verification evidence. Equations and citations are not decoration:
they make assumptions reviewable, help distinguish a published model from a code-specific
variant, and give future maintainers a route to investigate or improve the implementation.

For each substantial numerical algorithm or physical model, document the following where
applicable:

#. **Problem and assumptions:** governing equation or discrete problem, dimensionality,
   coordinate conventions, units, boundary/initial conditions, and validity limits.
#. **Mathematical definition:** define symbols and sign conventions; state the equation
   actually implemented, including limiting cases and relevant stability or convergence
   conditions. Do not copy an equation without checking it against the implementation.
#. **Source and provenance:** cite the original paper or book for the method, a current
   authoritative specification or maintainer resource when available, and the exact
   version or section used. Prefer DOI, publisher pages, official project manuals, or
   stable maintained URLs over search-result snippets.
#. **Implementation mapping:** identify the OneFLOW entry point and key functions, explain
   any differences from the cited formulation, and state whether a choice is a deliberate
   variant, a compatibility constraint, or an unresolved question.
#. **Verification:** link tests and representative cases to the properties they establish:
   analytic or manufactured solutions, conservation, limiting behavior, grid/time
   convergence, benchmark comparisons, or regression tolerances. A citation does not
   prove the code is correct, and a passing unit test does not establish model validity.
#. **Known limits:** record unsupported regimes, approximations, expected floating-point
   variation, and evidence that has not yet been collected.

Use references proportionally. A small utility function usually needs a contract and tests,
not a literature review. A flux scheme, turbulence closure, reconstruction method, limiter,
mesh-matching tolerance, or convergence algorithm may need a derivation or a clear citation
because its mathematical choices affect results. When a formula is standard but a particular
implementation choice is consequential, cite the standard formulation and explain the
choice locally.

Public open-source code and its theory/developer manuals can be valuable comparative
evidence, but they are not automatically authoritative. Record the project version or commit,
the relevant file or section, its license, and the exact behavior being compared. Treat
closed-source manuals as secondary references when they are lawfully accessible; cite their
title, publisher, edition/version, and relevant section without reproducing restricted text
or assuming undocumented implementation details. Do not copy code or substantial prose
from either open or closed sources without checking license and permission. Prefer an
independent derivation and OneFLOW-specific explanation over a line-by-line imitation of
another solver.

For CFD verification and validation, useful starting points include the
`NASA Turbulence Modeling Resource <https://tmbwg.github.io/turbmodels/>`_, which provides
model definitions and verification cases, and the
`NASA CFD Verification and Validation tutorial <https://www.grc.nasa.gov/www/wind/valid/tutorial/tutorial.html>`_.
Use a source appropriate to the actual model or algorithm; these resources are starting
points, not universal citations for every numerical method.

The intended maintenance loop is:

```text
theory/reference -> OneFLOW formulation -> implementation -> focused tests
        ^                                      |                |
        |                                      v                v
updated explanation <- discrepancy analysis <- case/regression evidence
```

When new evidence changes the understanding of an algorithm, update the algorithm page,
its tests or reference cases where warranted, and the source links together. Keep claims
labelled as verified, implemented-but-unverified, partial, or planned. This creates a
traceable improvement loop; it does not mean that an automated agent may silently change
scientific assumptions or accept its own output as proof. Changes to equations, model
variants, tolerances, or reference baselines require explicit review and reproducible
evidence.

Related design pages
--------------------

* :doc:`mesh_architecture` — mesh snapshots, canonical data, and serialization boundaries.
* :doc:`grid_topology` — topology invariants, boundary normalization, and derived connectivity.
* :doc:`workbench_architecture` — projects, workflows, run locations, artifacts, and process lifecycle.
* :doc:`info` — project scope and documentation entry points.
* `Documentation consolidation plan <https://github.com/eric2003/OneFLOW/blob/hexin/doc/plans/oneflow-documentation-consolidation-plan.md>`_ — document ownership and migration status.
* `Runtime refactoring roadmap <https://github.com/eric2003/OneFLOW/blob/hexin/doc/plans/runtime-refactoring-roadmap.md>`_ — staged refactoring goals and acceptance criteria.

Parameter schema and GUI editing
--------------------------------

See :doc:`case_parameter_design` for the sample-backed parameter classification and
for why a runtime database dump is not sufficient as an authoring model.

Directory separators and path semantics
---------------------------------------

Treat filesystem paths as paths, not as strings with a required final separator. Use
``std::filesystem::path`` composition (``root / relativePath``) when constructing child
paths, and use lexical normalization for comparisons that are intended to ignore redundant
separators or ``.`` components. Do not compare a path's raw string representation when the
contract is directory identity: a trailing separator may differ without changing the
represented directory.

OneFLOW currently has a compatibility wrinkle. ``Prj::ResolveCaseDir()`` deliberately
returns a directory string with a trailing separator because legacy call sites may rely
on that representation. Other path helpers already use ``std::filesystem::path`` and
normalize directory strings at selected boundaries. Therefore, the current contract is
not yet “all internal directory strings are separator-free.” Do not remove the trailing
separator globally until direct concatenations and legacy consumers have been audited and
migrated to path composition. The safe incremental rule is:

* preserve the existing public/legacy representation where compatibility requires it;
* compose paths with ``std::filesystem::path`` rather than string concatenation;
* normalize paths before semantic comparisons, while preserving filesystem roots;
* test both ordinary directories and root directories on supported platforms.

A directory-normalization helper must not strip a root such as ``/`` or ``C:\`` to an
invalid or semantically different string such as ``C:``. The root length is the boundary
for removing redundant trailing separators. This distinction matters on Windows, where a
drive-relative path (``C:folder``) and a drive-rooted path (``C:\folder``) do not mean
the same thing.

The runtime uses ``Prj::AreSameDirectoryPath()`` for the decision that selects the
default ``results/io_paths.txt`` location. This matters because ``ResolveCaseDir()``
retains a trailing separator while resolved output roots may not. Comparing their
``lexically_normal()`` paths directly can still treat these representations as unequal
and choose the wrong diagnostic location. The shared comparison helper applies the same
directory formatting policy used by ``WriteIOPathInfo()`` and rejects empty paths rather
than treating two missing values as a directory.

The decision can be summarized as:

.. code-block:: text

   configured mesh output root ----+
                                    +--> both roots equal case root? -- yes --> results/io_paths.txt
   configured flow output root -----+                            \-- no ---> io_paths.txt

Both mesh and flow output roots must identify the active case directory for the default
diagnostic to be written under ``results/``. An explicit ``ioPathInfoFileName`` still
overrides this default.

This is a staged compatibility policy, not a recommendation to keep trailing separators
forever. The long-term direction is to carry ``std::filesystem::path`` values internally
and convert to strings only at APIs that require them. Each migration should be justified
by call-site evidence and tests, not by changing formatting alone.


Control-file ingestion and its remaining module boundary
---------------------------------------------------------

The production case-specific entry point is
``ReadControlInfo(caseDir)`` in ``codes/io/src/ControlFileReader.cpp``. Its
responsibility is orchestration rather than numerical interpretation:

#. Read the case's ``script/control.txt`` manifest.
#. Resolve each listed script under the case directory and parse the files in order.
#. Commit the resulting typed entries into the global ``DataPara`` store.
#. Coordinate process status and broadcast the committed parameter database.
#. Write the merged database snapshot to ``log/database.log``.

The configuration pipeline already contains a useful internal separation:
``ConfigLoader`` parses legacy text into a ``ConfigDocument``, while
``ConfigDatabaseAdapter`` commits entries using their declared type names.
The database then owns the typed values. A test of the real control reader must
exercise these production stages rather than instantiate a test-only parser.

The source-level split is now implemented. Reader orchestration is in
``ControlFileReader.cpp``; script parsing and ConfigDocument commit are in
``ControlScriptParser.cpp``; runtime array-size lookup is in
``ControlParameterResolver.cpp``; database broadcast, compression/
decompression, and diagnostic dumping are in ``ControlDatabaseIO.cpp``;
legacy parameter inspection helpers remain in ``ParaFile.cpp``. Public
declarations in ``ParaFile.h`` are unchanged. The repository's ``ConstructSolutionDirTree()``
source collector recursively discovers C++ sources, so these files are expected
to enter the production source list without a per-file CMake registration.
This expectation still requires an exact-revision build.

The split reduces responsibility overlap, but it does not yet make the real
reader independently linkable: ``ReadMultiScriptFiles()`` calls
``ReadOneFLOWScriptFile()`, which remains in ``ParaFile.cpp``. The next
dependency review must determine whether to extract this parser adapter and its
minimal dependencies into a third production unit or expose an existing narrow
target. Do not duplicate parsing in a test.

The component integration test must initialize the real serial parallel backend
through ``SetUpParallelEnvironment()``, use a temporary multi-file control
manifest, verify typed database values and the case-local database dump, and
clear global database state afterwards. Sequential-case isolation and setup
failure cleanup remain separate lifecycle acceptance criteria.

The system-level benefit is a clearer ingestion pipeline and a first seam for
replacing legacy syntax adapters without changing consumers. It does not alter
any numerical algorithm or physical model. R1 remains partial until production
reader integration, effective case path application, sequential-case isolation,
and setup-failure cleanup have reproducible test evidence.

For the broader target and staged evolution model, see
:doc:`evolutionary_computation_architecture`.


The resolver boundary has a focused test target,
``control_parameter_resolver_test``, which exercises literal lengths,
named lengths through the real DataBase, and missing-name failure. The
``ConfigLoader_test`` suite separately exercises parser behavior with an
injected resolver. These two test levels deliberately distinguish parsing
semantics from runtime lookup. Neither replaces the pending production
multi-file reader integration test.


Case teardown and cleanup failure semantics
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

SimuContext::TeardownCase releases case-owned resources in dependency order:
accelerator-backed state and shared backend buffers first, solver and physics
state next, then zone/field/database state, and finally case log and
project-path bindings. The process-level accelerator and parallel runtimes
remain alive until the process owner calls FinalizeEnvironment().

Each release/reset is an independent cleanup step. The production teardown
uses SimuCaseLifecycleDetail::ExecuteAllCleanupSteps() so an exception in
one step does not prevent later steps from being attempted. The helper records
the first cleanup exception, executes all remaining steps in order, and then
rethrows that first exception. Context-owned task metadata is reset as the
final step even if an earlier release failed.

When teardown runs from SimuImp::RunCase() while a task/setup exception is
already active, ExecuteWithCleanup() preserves the original execution
exception; teardown failures do not replace it. This policy guarantees that
all cleanup attempts are made, not that every release succeeds. A successful
case whose teardown fails still reports the first teardown exception. Focused
unit tests pin ordering, continued cleanup after multiple failures, and
first-failure preservation. These tests do not replace an integration test of
real solver/device resource release.


Case environment readiness
--------------------------

``SimuContext::IsEnvironmentReady()`` reports successful completion of the
entire case bootstrap, not merely that setup has started. The flag is cleared
at the beginning of ``SetupCaseEnvironment()`` and set only after the case
root/log binding, production control ingestion, IO-setting validation, and
IO-path application have all completed. If any step throws, the owning
``SimuImp::RunCase()`` or compatibility rollback path still attempts case
teardown, while the context remains not-ready. This distinction prevents
callers from observing a partially initialized case as ready. It does not
replace the outstanding integration test of the full production bootstrap.

Production case-bootstrap composition
-------------------------------------

The ordered case bootstrap is factored into
``SimuContextDetail::InitializeCaseConfiguration(caseDir)`` in
``codes/main/src/SimuContextBootstrap.cpp``, declared by
``SimuContextBootstrap.h``. This is the production function called by
``SimuContext::SetupCaseEnvironment()``; it binds the case project root and
log, invokes the real ``ReadControlInfo(caseDir)``, validates the typed IO
settings snapshot, and applies the effective generic and mesh/flow roots.
The serial integration target links this same production translation unit
and exercises two successive case configurations, including configured
roots on the first case and default-root rebinding on the second. A scoped
test guard clears case-owned log, project-path, and parameter state even when
a fatal assertion exits the test early. The test still resets state directly
between the two bootstrap calls; it does not invoke production
``TeardownCase()``, the full ``SimuImp::RunCase()``, MPI collectives, or
accelerator resource release.

This is meaningful component-integration evidence: the production bootstrap
composition, real control reader, typed database commit, IO-setting
application, and sequential root rebinding are exercised together. R1
remains partial because the public readiness transition in
``SetupCaseEnvironment()`` and cleanup wiring through production
``RunCase() -> PostProcess() -> TeardownCase()`` are not yet exercised
end-to-end.

Process runtime ownership and Context moves
-------------------------------------------

SimuContext is move-constructible because a context may be transferred to
another owner. Its process-runtime marker must move with it: copying the
Active phase would leave both the source and destination claiming to own
the same process lifecycle. ProcessRuntimeState is therefore non-copyable
and has a transfer move constructor. The moved-from marker enters a terminal
MovedFrom phase, so it cannot be started again accidentally. Move
assignment is disabled for both the marker and SimuContext because
overwriting a destination could abandon an already active process-runtime
owner. This is an ownership rule, not a resource-finalization operation; the
destination remains responsible for the eventual explicit process teardown.


Moved-from SimuContext contract
-------------------------------

The explicit SimuContext move constructor transfers its arguments, case
identity, task selection, readiness marker, process lifecycle marker,
accelerator-state registry, and expanded solver names. It then resets the
source context's observable scalar state and clears its remaining metadata.
In particular, a moved-from context must not continue to report that its case
environment is ready or its task is resolved. Move assignment remains
disabled because replacing an existing context could abandon owned runtime
state. This move operation transfers ownership; it does not perform resource
teardown or process finalization.


Accelerator backend-state ownership
------------------------------------

``EulerDomainStateRegistry`` owns backend state objects through unique ownership.
Its container is held behind a ``std::unique_ptr<StateMap>``, allocated lazily on
the first insertion. This is intentional: the registry's move operation transfers
one pointer and therefore does not depend on the standard library's implementation-
specific ``noexcept`` traits for moving or default-constructing
``std::unordered_map``. The backend state objects themselves are not moved.
After a move, the source registry is empty and remains valid; a later insertion
allocates a fresh map.

This is a small indirection with a clear lifecycle benefit. The registry is a member
of the simulation context, whose move contract must not throw after process-lifecycle
ownership has begun transferring. The invariant is tested at compile time and by
a runtime ownership-transfer test. The trade-off is one allocation on first use and
one pointer indirection for registry operations; an empty registry performs no
allocation.
