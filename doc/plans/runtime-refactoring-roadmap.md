# OneFLOW Runtime Refactoring Roadmap

Status: In progress
Branch baseline: `hexin`
Last reviewed: 2026-10-30

This is a living engineering roadmap, not a promise that every proposed boundary already
exists. It records the direction and acceptance criteria for incremental runtime refactoring.
The source code, tests, and CI results for the exact commit remain the authority for
implementation status.

## 1. Working method

For substantial stages, start from the reusable [refactoring stage record template](templates/refactoring-stage-record.md). Keep the record focused on durable engineering knowledge; the [AI project context](../AI_CONTEXT.md) is the short current-state index, not a duplicate of every stage record.

For each refactoring slice, keep a short record of:

1. **Observed problem** - a concrete behavior, maintenance cost, dependency, or unclear contract.
2. **Current path** - entry point, direct callers, state ownership, and lifecycle.
3. **Intended contract** - what must remain true and what behavior should change.
4. **Smallest useful seam** - a function, value object, adapter, or interface only when it makes the contract easier to test or maintain.
5. **Evidence** - focused tests, relevant integration/regression tests, and CI links tied to a commit.
6. **Outcome** - completed, partially completed, blocked, or deferred; include follow-up gaps.

Do not record every command or private brainstorming step. Keep durable reasoning, decisions,
and verification evidence. A refactor is not complete merely because code moved into a new
class or all tests unrelated to the changed contract passed.

## 2. Cross-cutting constraints

- Work from the latest `hexin` HEAD; do not target `master`.
- Preserve behavior unless a behavior change is explicitly intended and documented.
- Follow actual call sites and ownership before changing interfaces.
- Keep parsing, validation, and domain interpretation distinct where the current code supports that separation.
- Prefer deterministic tests at the narrowest meaningful boundary, followed by integration and regression checks.
- Keep C++ source files ASCII-only.
- Update the curated design page when a contract or responsibility changes. Add an ADR only for decisions with durable consequences or meaningful alternatives.
- Never describe a planned feature or unrun test as implemented or verified.
- Keep `doc/AI_CONTEXT.md` as the fast handoff entry for AI-assisted work; update it at material milestones with a verified base commit, current focus, known gaps, evidence, and next actions.

## 3. Runtime refactoring sequence

The sequence is provisional and should be adjusted when code inspection or test evidence changes
the risk assessment. Finish one independently verifiable slice before expanding its scope.

### R0. Establish a source-backed map

**Purpose:** Make entry points, execution phases, configuration layers, and known limitations
visible before selecting local refactors.

**Current evidence:** `doc/source/runtime_architecture.rst` maps
`OneFlow.cpp -> Simulation -> SimuImp -> SimuContext -> TaskRegistry -> ISimuTask`
and distinguishes process and case lifetimes.

**Exit criteria:**
- The map matches the current `hexin` source tree.
- Current behavior and proposed boundaries are explicitly distinguished.
- The next slice is selected from a concrete call path, not from a desire to introduce a pattern.

**Status:** Initial map documented; continue validating it as code changes.

### R1. Make case configuration and I/O root resolution understandable

**Observed area:** `SimuContext::SetupCaseEnvironment()` in
`codes/main/src/SimuContextEnv.cpp` reads optional root settings, applies project and
generic I/O roots, applies mesh/flow-specific roots, falls back to `sourceGridCaseDir`
for mesh input, and writes path information.

**Questions to verify before changing code:**
- Does parameter lookup distinguish a missing key from a present key with an invalid type or value?
- What is the exact base for each relative path, and in what order do generic and domain-specific roots apply?
- Which empty and `default` values preserve compatibility?
- Which globals or filesystem side effects make setup difficult to test?
- What cleanup is required if setup fails partway through?

**Intended approach:**
1. Trace `GetDataValue` and `Prj` path setters plus direct callers and existing tests.
2. Write behavior tests for the verified precedence and fallback rules before extraction.
3. If the code supports it, isolate deterministic root resolution from global mutation and file output.
4. Keep application of resolved paths and writing `io_paths.txt` at the environment boundary.
5. Update the runtime architecture page with the verified contract and record test/CI evidence.

**Exit criteria:**
- Existing case-relative behavior and `sourceGridCaseDir` compatibility are covered by tests.
- Missing settings are not accidentally conflated with invalid settings unless that is the established contract.
- Resolution rules can be tested without initializing the full solver where practical.
- Relevant integration tests and Windows/Linux CI are reported against the changed commit.

**Status:** Partial. Root-setting lookup and legacy fallback selection are isolated into a tested value-object adapter. The real `SetupCaseEnvironment()` / `ReadControlInfo(caseDir)` integration path, effective `Prj` application, diagnostic output, and sequential-case behavior remain unverified.

### R2. Trace task selection and parameter semantics

**Purpose:** Establish which layer parses control syntax, which layer validates parameter values,
and which layer applies task-specific meaning.

**Candidate path:** `ConfigLoader`, the runtime control database, and
`SimuContext::ResolveTaskFromControl()`, followed by the registered `ISimuTask`.

**Exit criteria:**
- A representative task path is documented end to end.
- Required versus optional parameters and default behavior are supported by code/tests.
- Parsing changes do not silently change task semantics.
- Any extracted adapter has a concrete consumer and focused tests.

**Status:** Partial. `ConfigLoader_test` covers syntax-to-document behavior; runtime parameter validation and task-specific meaning remain open.

### R3. Trace one grid task and one solve task

**Purpose:** Compare the intended task/context boundary with actual dependencies on free functions,
global state, solver registries, and runtime resources.

**Exit criteria:**
- One grid-generation path and one solve path have concise call-path notes.
- Ownership and teardown responsibilities are identified.
- A seam is extracted only when it enables a meaningful test or removes a real dependency.
- Numerical behavior is protected by applicable grid/CFD regression evidence.

**Status:** Planned; scope depends on R2 findings.

### R4. Align Workbench execution contracts with runtime behavior

**Purpose:** Verify project/case selection, run location, workflow dependencies, artifact transfer,
history, and process lifecycle against implementation and tests.

**Exit criteria:**
- Independent multi-case execution is not confused with dependency workflows.
- Source case inputs and staged run directories follow explicit contracts.
- Run records reflect actual process outcomes and locations.
- Background-task behavior is documented as implemented only after code and tests prove it.

**Status:** Planned; use `doc/source/workbench_architecture.rst` as the curated overview.

### R5. Consolidate stable contracts into user and developer guidance

**Purpose:** Turn verified architecture into practical documentation without duplicating every
engineering note.

**Exit criteria:**
- User instructions match tested command behavior and paths.
- Developer guidance links to the architecture map, roadmap, relevant runbooks, and evidence.
- Historical documents remain clearly identified until their content and inbound links are reviewed.
- Documentation build and relevant CI evidence are tied to the current revision.

**Status:** Ongoing alongside earlier stages.

## 4. Per-slice progress record

### 2026-10-30 - Read all case I/O settings before path mutation

- Based on inspected `hexin` HEAD `df567dba6f77abba725f05ce72cb32499046a653`; user reports 967/967 tests passed on the pre-change baseline.
- Expanded the root-only `CaseIORootSettings` snapshot into `CaseIOSettings`, including `ioPathInfoFileName`. Production setup now reads and type-checks that setting before applying configured roots to global `Prj` state.
- Added the diagnostic filename to the invalid-type contract test and full snapshot test. This verifies the configuration snapshot contract, not the complete production setup path.
- Updated the runtime guide with a setup-order diagram and the limitation that the snapshot separates reads from mutations but does not make path application or filesystem output transactional.
- Added a ranked whole-system improvement sequence and a staged source-comment-to-Sphinx documentation strategy to the architecture review and runtime guide.
- Validation: source changes checked for expected symbols and ASCII-only C++ content. No C++ build, test run, documentation build, or exact-commit CI result is available in this environment.
- System-level assessment: all applicable setup-time string settings are selected before root mutation, reducing partial configuration application on a type error. R1 remains Partial because real control-reader integration, effective path binding, diagnostics, sequential-case isolation, and failure cleanup are not covered end-to-end.


Append a short entry when a slice materially advances. Keep it factual and link the exact commit.

### 2026-10-30 - Extend the R1 root-setting type contract

- Based on inspected `hexin` HEAD `5ae294d5956e8811762d96b402e93f3a77fece5b`; the user reports 966/966 tests passed on the pre-change baseline.
- Extended `simu_context_parameter_test.cpp` so every root-setting key, including the conditional legacy `sourceGridCaseDir`, is checked for rejection when its declared type is not a string.
- This is a narrow schema-contract test. It does not claim to exercise `ReadControlInfo(caseDir)`, effective `Prj` root application, path diagnostics, sequential-case reset, or failure cleanup.
- Validation: test source and current call path inspected; no C++ build or test run in this environment. Exact-revision CI has not been checked.
- System-level assessment: this improves the boundary between configuration interpretation and global path mutation, but does not change R1's Partial status. Next action remains a narrow production-owned setup orchestration seam that preserves the real control reader.

### 2026-10-26 - Initial roadmap

- Added a staged runtime refactoring sequence and acceptance criteria.
- Selected case I/O root resolution as the first code investigation because setup currently combines configuration reads, root precedence, global path application, and path-info output.
- No C++ behavior change is included in this roadmap entry.
- Validation: documentation and runtime CI were still in progress at the recorded baseline; this entry itself does not claim a passing build or test suite.

## 5. Review checklist

Before calling a slice complete, ask:

- Is the original problem stated in terms of behavior or a concrete maintenance burden?
- Have direct callers, state ownership, path semantics, and failure cleanup been checked?
- Do tests distinguish the intended contract from a plausible wrong implementation?
- Is the code change smaller than the problem boundary and free of unused abstractions?
- Are compatibility behavior and intentional changes explicit?
- Are test commands, results, commit SHA, and CI links accurate?
- Do the design page, roadmap status, and user-facing instructions agree with the code?


### 2026-10-26 - R1 call-path inspection

Verified from the current source:

- `SimuImp::InitSimu()` calls `SimuContext::SetupCaseEnvironment()`; the latter currently binds case state, reads optional root settings, applies the generic roots, applies domain-specific roots, and writes path diagnostics.
- `Prj::SetPrjBaseDir()` resets project and domain-specific roots to prevent settings leaking between cases.
- Generic and domain roots are applied in two steps. Empty or `default` domain roots inherit generic roots; empty, `default`, or `case` generic roots resolve to the active case. The special `project` root requires a configured project root. Relative roots other than these tokens are resolved relative to the active case; absolute roots are used directly.
- The mesh input fallback from `sourceGridCaseDir` is applied only when `meshInputRoot` is empty or `default`. An explicit mesh input root takes precedence.
- `GetDataValue(name)` throws when a parameter name is absent. However, the `DataObject` typed lookup uses a static cast to the requested element type; it does not provide a safe runtime type check. Therefore, broad `catch (const std::exception&)` blocks reliably treat a missing key as optional, but they must not be described as robustly handling wrong-typed values. A wrong type is not guaranteed to throw and may be unsafe.
- Existing `tests/main/prj_cmdline_test.cpp` covers many `Prj` root and file-resolution rules, while `tests/main/simu_context_test.cpp` primarily tests context state and task plumbing. The inspected tests do not directly exercise the full `SetupCaseEnvironment()` parameter precedence and legacy fallback path.

**Updated next action:** do not extract a resolver yet. First add or identify a test seam that can verify root precedence and `sourceGridCaseDir` fallback without booting the complete runtime. In parallel, inspect how parameter entries are constructed and typed so the test uses valid database values. Only after this contract is covered should the code boundary be changed.

Validation note: this entry records source/test inspection only. No C++ source or test was changed and no local test suite was run for this entry.


### 2026-10-26 - R1 first implementation slice

Commit: `d8568fe91722af5d78d8083d6bae23d9446d952d`

- Replaced broad exception swallowing for the optional path-related string parameters with a small local reader that distinguishes a missing entry from a present entry whose declared database type is not `HX_STRING`.
- Missing optional entries still resolve to an empty string, preserving the current fallback path. Present non-string entries now raise `std::invalid_argument` instead of being reinterpreted through an unchecked `static_cast`.
- Applied the same rule to `sourceGridCaseDir` and `ioPathInfoFileName`; the explicit `meshInputRoot` precedence and legacy source-grid fallback condition are unchanged.
- The helper remains local to environment bootstrap. No general database API or new public abstraction was introduced.
- Validation: source is ASCII-only by inspection. CI for this commit must still confirm compilation and regressions. A focused automated test for this specific invalid-type setup path has not yet been added, so this is a first implementation slice, not a completed R1 stage.


### 2026-10-26 - AI-assisted development handoff

- Added `doc/AI_CONTEXT.md` as the fast entry point for restoring project state across AI sessions and developer handoffs.
- Linked the handoff from the repository README, documentation README, and Sphinx developer resources so it is discoverable without reading the full repository.
- The page summarizes branch constraints, the current R1 focus, source-backed path behavior, the last implementation slice, explicit verification gaps, and prioritized next actions. It requires rechecking the current `hexin` HEAD and CI before relying on the snapshot.
- This is a documentation/process improvement only; it does not claim any new code behavior or test result.


### 2026-10-26 - Reusable stage record and include cleanup

- Added `doc/plans/templates/refactoring-stage-record.md` so substantial refactoring stages use a consistent record for problem, call path, intended contract, minimal boundary, changes, evidence, open gaps, and next action.
- Linked the template from this roadmap and `doc/AI_CONTEXT.md`; the handoff page remains a concise overall state index rather than a duplicate of every stage.
- Added an explicit `<stdexcept>` include in `codes/main/src/SimuContextEnv.cpp` for the existing `std::invalid_argument` use (commit `ec9575216e9b6db14a400446c1373e22405c7169`).
- Validation: this was a source/include correctness cleanup and documentation change. No local build or tests were run in this step; the invalid-type behavior still lacks a focused automated test. Continue R1 test-seam investigation and verify CI for the exact current branch revision.


### 2026-10-26 - Whole-system architecture review mechanism

- Added `doc/source/architecture_review.rst` as a handbook chapter and navigation entry.
- The chapter provides a reusable CFD/CAE concern map, calls out OneFLOW cross-cutting risks (global runtime state, configuration validation, mesh derived-data invalidation, workflow artifact provenance, lifecycle cleanup, and evidence discipline), and includes a repeatable prompt for AI-assisted engineering.
- Linked the chapter from `doc/AI_CONTEXT.md` and `doc/README.md`. The intended workflow is to zoom into a concrete call path, implement and test a small change, then zoom back out at meaningful milestones and update priorities.
- This documentation change does not resolve the R1 runtime parameter test gap. R1 remains partial / implemented but unverified until the relevant contract has focused test and exact-revision validation evidence.
- Validation: no Sphinx build, C++ build, unit test, or CI result is claimed by this documentation update.


### 2026-10-26 - Extract and test optional string parameter lookup

- Added `codes/main/include/SimuContextParameter.h` with a small helper that accepts `const DataPara&`. This keeps the parameter type contract directly testable without initializing the full simulation runtime or adding a public database API.
- Updated `codes/main/src/SimuContextEnv.cpp` to use the helper while retaining the existing global-database lookup and call-site behavior.
- Added `tests/main/simu_context_parameter_test.cpp` and registered `simu_context_parameter_test` in `tests/main/CMakeLists.txt`. The tests cover a missing optional entry returning an empty string, reading a valid string, and rejecting a non-string entry with the parameter name and expected type in the message.
- Implementation commits: [helper](https://github.com/eric2003/OneFLOW/commit/fd3889362408586e379d45d08a170abf2077efb9), [production integration](https://github.com/eric2003/OneFLOW/commit/6d72e31f93f55b00b622e122fa75eb689832694f), [tests](https://github.com/eric2003/OneFLOW/commit/1e438f84f600c689666f831ff0837ce3d5936e13), [CMake registration](https://github.com/eric2003/OneFLOW/commit/a13a890bb31e2d84434fbe4297cf76aebc24c865).
- Validation: source and CMake files were inspected after the changes. The exact-revision GitHub combined status query returned no checks, and no local C++ build or test execution was available in this step. Do not interpret the new tests as passing until `simu_context_parameter_test` is built and run.
- Remaining R1 work: exercise path precedence, `sourceGridCaseDir` fallback, and per-case path reset through a realistic setup-level test. The helper-level contract is implemented but unverified; R1 remains partial.


### 2026-10-26 - R1 evidence update and runtime guide correction

- The user reports `100% tests passed out of 697` on the latest `hexin` code baseline. This is user-reported local evidence; no independent build/test was run in this assistant session, and no exact-revision GitHub CI pass is claimed.
- Updated `doc/source/runtime_architecture.rst` to remove the stale statement that optional path lookups still use broad exception handlers. It now documents the actual `SimuContextDetail::ReadOptionalStringParameter()` type check and distinguishes helper-level coverage from full setup integration.
- Updated `doc/AI_CONTEXT.md` with the reported test count and the remaining acceptance gap.
- System-level assessment: the type contract for optional path strings has a narrow unit seam, and existing `Prj` tests already cover generic/domain root behavior, inheritance, diagnostics, and root reset. The remaining concrete gap is testing how `SetupCaseEnvironment()` applies the control database values, especially the explicit `meshInputRoot` versus `sourceGridCaseDir` fallback. Do not create another path resolver until that behavior is exercised.
- Status: R1 remains partial. Next step is to inspect `ReadControlInfo()` and global database initialization to determine whether a realistic temporary-control-file integration test can be added without reproducing the whole solver bootstrap.


### 2026-10-26 - Isolate and test the source-grid fallback decision

- Extracted the condition that decides whether `sourceGridCaseDir` should be read into `SimuContextDetail::ShouldUseLegacySourceGridCaseDir()` in `SimuContextParameter.h`. Production setup now calls this helper, while the source-grid value remains lazily read only when the mesh input root is empty or `default`; this preserves the existing behavior that an explicit mesh root does not trigger a read of the legacy fallback parameter.
- Added focused tests: empty and `default` mesh roots permit fallback; explicit, `case`, and ordinary relative mesh roots do not select the legacy fallback. Renamed the predicate to `ShouldUseLegacySourceGridCaseDir()` because it expresses the policy decision rather than the mechanics of reading a parameter.
- Implementation commits: [fallback policy](https://github.com/eric2003/OneFLOW/commit/f4731c5b273cc9774e02b40a97ba49e66411c6f3), [production call site](https://github.com/eric2003/OneFLOW/commit/afd6cd2d3eaa1fed8b35f1705f36f2287fd62ce5), [tests](https://github.com/eric2003/OneFLOW/commit/897112bde18f9f7ab77b063605bd56f3b0e0f749).
- Validation status: the user-reported 697/697 pass predates this newest code/test change. The added tests and latest documentation commits have not yet been run in a fresh local suite in this session; do not claim they passed.
- Updated the runtime architecture guide to document the checked optional-string helper and corrected a malformed literal line break in its related-links list.
- Remaining gap: the new test validates the branch policy, not the complete control-file-to-`SetupCaseEnvironment()` path. `ReadControlInfo(caseDir)` performs parallel coordination, control-parameter broadcast, and database dump, so any integration test should reuse existing runtime bootstrap rather than mock only the branch or reproduce the entire solver initialization.


### 2026-10-26 - Clarify legacy mesh fallback policy naming

- Based on the latest inspected `hexin` baseline, renamed the predicate to `ShouldUseLegacySourceGridCaseDir()`. The old name described an I/O operation (read) rather than the policy (use fallback).
- Added negative contract tests for the special `case` token and a normal relative mesh root. The existing positive tests for empty and `default`, plus the explicit-root test, remain.
- Commits: [policy name](https://github.com/eric2003/OneFLOW/commit/d4e858abe3aa1a37fa8b6a11b8c6cd50ebf13624), [production call site](https://github.com/eric2003/OneFLOW/commit/fa73d7dbb04e99537a82cfd2aaf330837815607f), [contract tests](https://github.com/eric2003/OneFLOW/commit/8cf20d872d844b0fa68d99c5a59661be95ef0192).
- Validation: source and test files were re-fetched from GitHub and checked for ASCII-only content. No build, test execution, or successful CI status is available for this revision. The user-reported 700/700 pass predates these changes.
- Architecture assessment: the helper-level test is cheap because it isolates a deterministic policy; it still does not prove the complete `SetupCaseEnvironment()` flow. Keep R1 partial until an integration seam is found and verified.

### 2026-10-26 - Test the effective legacy mesh input fallback value

- The previous tests covered whether fallback lookup is allowed, but did not test the effective value chosen when the legacy parameter is empty, `default`, or a valid source case. This left the main compatibility outcome implicit in the production `if` block.
- Added `ApplyLegacySourceGridCaseDir()` as a small deterministic policy function and used it at the production call site. Parameter lookup remains lazy: `sourceGridCaseDir` is still read only when `meshInputRoot` is empty or `default`.
- Added tests for empty/default mesh roots, empty/default legacy values, explicit relative roots, and the `case` token. These tests describe the value contract without requiring the MPI/control-file bootstrap.
- Commits: [policy helper](https://github.com/eric2003/OneFLOW/commit/d643aed3f964d9e9ea48947a259dd6e0820bd732), [production call site](https://github.com/eric2003/OneFLOW/commit/e38adea4fa70e51299cb14da3164209fcd89fd4c), [contract tests](https://github.com/eric2003/OneFLOW/commit/8299d15eb1241f60f949dc1b081ef78cf32f905f).
- Validation: source and tests are committed, but no local build/test run or successful exact-revision CI result is available in this session. Do not count these tests as passed until they run.
- Architecture assessment: the first helper tested only a branch predicate; adding the effective-value test exposed a better contract boundary. The setup integration gap remains because `ReadControlInfo(caseDir)` performs script loading, process coordination, parameter broadcast, and database dump. Keep R1 partial.
### 2026-10-26 - Extract and register case cleanup contract

- Inspected `hexin` HEAD `b30a86ed`. `SimuImp::RunCase()` now delegates the execute/cleanup exception policy to the production-owned `SimuCaseLifecycleDetail::ExecuteWithCleanup()` boundary.
- The contract is: cleanup runs once after successful execution; cleanup also runs after execution failure; an execution failure is preserved if cleanup also throws; a cleanup failure is propagated when execution succeeds.
- The focused test target is registered in `tests/main/CMakeLists.txt`. This session strengthened the success-path assertion to verify ordering, not only callback counts.
- Validation: the existing Windows build directory lists 618 CTest entries, with some stale or unavailable targets represented as `NOT_BUILT`. CMake began regeneration for `simu_case_lifecycle_test` but did not produce a completed build or executable in this session. No test pass or 700/700 result is claimed for `b30a86ed`.
- Architecture assessment: the cleanup policy can be tested without linking the solver runtime, which confirms a small lifecycle seam is useful. It does not validate full `SimuContext::TeardownCase()` behavior or the real control-file setup failure path. R1 remains partial.
- Next: complete focused Windows target validation, then revisit a production-linked temporary-case integration boundary for `SetupCaseEnvironment()` without duplicating executable source configuration.

### 2026-10-26 - Make multi-case process cleanup exception-safe

- In `Simulation::RunImpl()`, process setup previously happened before the `try` block. An exception from accelerator initialization therefore bypassed `FinalizeEnvironment()`. The case loop did finalize on success and failure, but duplicated cleanup handling.
- Added `RunCasesSequentially()` as the production policy boundary and wrapped process setup plus the ordered case loop with `ExecuteWithCleanup()`. The cleanup callback now runs after success, case failure, or process setup failure; the existing helper preserves the original execution exception if finalization also throws.
- Added tests using the user's six case paths to assert exact ordering, plus a failure test proving later cases are skipped and process finalization runs.
- Validation: the user reports 712/712 local tests passed on the prior baseline. The current Windows build folder lists 618 tests and includes `NOT_BUILT` entries. All seven lifecycle tests passed after direct compilation with the installed MSVC toolchain, and `Simulation.cpp` compiled separately with the project include paths and Windows definitions. The generated CMake/MSBuild target has not rebuilt the latest seven-test file, and the full executable link remains unverified. `Checking File Globs` is the generated CMake `CONFIGURE_DEPENDS` check; the current snapshot had no active CMake or VS build process, so concurrent CMake GUI and Visual Studio use of the same build folder remains a plausible but unconfirmed source of contention.
- Architecture assessment: process and case lifetimes are now both represented at the orchestration boundary. The tests exercise ordering and cleanup policy without booting MPI or accelerator runtime; a full six-case solver run remains separate integration evidence.

### 2026-10-26 - Close standalone single-case process lifecycle

- `SimuImp::Run()` initialized process runtime and executed `RunCase()`, but did not call `FinalizeEnvironment()` on success or failure. The multi-case entry point already owned process finalization, leaving the public single-case compatibility entry point with a different lifetime contract.
- Added `ExecuteWithProcessLifecycle()` to share setup / execution / finalization policy. `SimuImp::Run()` now uses it; multi-case orchestration composes through the same helper.
- Added a focused test that verifies finalization occurs after a failed single-case execution and the original exception is preserved.
- Validation: eleven focused lifecycle and process-state tests passed after direct MSVC compilation. `SimuImp.cpp` and `SimuContextEnv.cpp` compiled separately with the project include paths, Windows definitions, C++20, and MPI headers. The user reports 714/714 tests on the prior baseline; the full suite and executable link were not rerun for this slice.
- Architecture opinion: process lifetime and case lifetime should remain separate owners. The single-case path now owns both in sequence, while `RunCase()` remains reusable by the multi-case owner. A possible future concern is that per-case contexts still default their `Rank()` / `Size()` to 0/1; no current production call sites use these accessors, so propagate process values when a real consumer or test seam requires it.

- User requirement: MPI process initialization/finalization is one-shot, while the same case may be run more than once. The multi-case loop preserves every command-line entry, including identical case paths; each execution gets fresh case state, and the shared process context owns MPI lifetime. A single `RunCase()` call does not initialize or finalize MPI. Eleven focused tests validate repeated paths, setup/run/finalize ordering, exception cleanup, and the process lifecycle state.

### 2026-10-26 - Make compatibility setup and teardown failure-safe

- The compatibility ``SimuContext::SetupEnvironment()`` combined process startup and case setup without rollback if case setup failed. ``TeardownEnvironment()`` skipped process finalization if case teardown threw.
- Added ``ExecuteWithRollback()`` for setup-only operations: successful setup leaves resources active, while failed setup attempts rollback and preserves the original exception. Process setup uses it after claiming the one-shot lifecycle; combined setup uses it around case setup. Teardown composes ``ExecuteWithCleanup()`` so process finalization is attempted even when case teardown fails.
- Added tests proving rollback does not run on success and runs on failure without replacing the setup exception.
- Validation: all 13 focused lifecycle/state tests pass after direct MSVC compilation. ``SimuContextEnv.cpp`` compiles separately with MSVC, C++20, project includes, and MPI headers. The user reports 719/719 tests on baseline ``23e11c91``; the complete suite and executable link were not rerun for this slice.
- Architecture assessment: this improves partial-initialization and cleanup behavior at an existing compatibility boundary without adding a new runtime abstraction. Next priority remains a production-linked setup fixture and then separating process ownership from the case-oriented ``SimuContext`` type, guided by actual dependencies rather than a broad rewrite.

### 2026-10-26 - Add a minimal reusable CLI parsing core

- The existing parser lived behind ``Prj::ParseCmdLineArgs()`` and exposed no reusable argument-schema contract. OneFLOW currently needs a positional mode and one-or-more case directories, including repeated identical paths.
- Added header-only ``CliParser`` with declarative positional cardinality, named parsed values, schema validation, readable validation errors, and generated usage. The legacy ``Prj`` parser now adapts its existing argv vector to that schema; its public result type and mode semantics remain unchanged.
- Scope is deliberately limited to current use: no named switches, subcommands, config-file merge, callbacks, or automatic typed conversions. Add those only when a concrete OneFLOW command requires them.
- Added six parser-core tests and two adapter regressions for usage on missing values and preservation of duplicate case paths. The existing ``prj_cmdline_test`` target passed 46/46 when built directly with MSVC. The user reports 721/721 on baseline ``cdc51377``; the complete suite has not been rerun. No CMakeLists file changed, but the root source-tree glob will detect the new header and may trigger one CMake regeneration/reload in Visual Studio.
- Migration plan: retain ``Prj::ParseCmdLineArgs()`` while call sites move to a dedicated run-options value API; then assess mode validation and named options separately. Keep parsing independent of global ``Prj`` path state.

### 2026-10-26 - Migrate production callers to canonical run options

- The first CLI slice added a reusable parser but left production consumers reading ``CmdLineOptions``, which duplicated the first case in ``caseDir`` and the full list in ``caseDirs``.
- Added ``SimulationCommandLineOptions`` with a single ``caseDirectories`` collection and a free ``ParseSimulationCommandLineOptions()`` API. ``Simulation``, ``SimuContext``, and ``Prj::ProcessCmdLineArgs()`` now use it. ``Prj::ParseCmdLineArgs()`` remains a compatibility adapter for existing callers and tests.
- Existing argument count, mode interpretation, case order, and duplicate path behavior are unchanged. Top-level selector/run mode classification remains in ``Simulation::RunImpl()`` for a separate future contract change.
- Validation: the expanded ``prj_cmdline_test`` target passed 48/48 with direct MSVC compilation; ``Simulation.cpp`` and ``SimuContextEnv.cpp`` also compiled separately. The user-reported 721/721 preceded the original CLI parser addition; no full-suite result is claimed for this revision.
- Architecture assessment: this removes the duplicated case-directory state from production consumers while preserving the compatibility boundary. Next inspect invocation-mode branching and define a testable application command model before adding flags or changing syntax.

### 2026-10-26 - Make top-level invocation selection testable

- Starting from ``hexin`` HEAD ``9af8dc3a343fb6fa1192b4a04710b847d1d17be0``. Extracted the existing ``Simulation::RunImpl()`` argv-count branches into ``ClassifySimulationInvocation()`` and a named ``SimulationInvocationKind`` enum. The caller still owns execution; parsing full-run values remains a separate responsibility.
- Preserved existing behavior: one argv entry selects lightweight test mode, two select the current wrong-argument-count error, and three or more enter full simulation. Added one focused test for each path.
- Validation: direct MSVC compilation and execution of ``prj_cmdline_test`` passed 51/51; ``Simulation.cpp`` compiled separately with project include paths and MPI headers. No full-suite result is claimed for this slice; the user-reported 721/721 predates CLI parser changes.
- Architecture assessment: mode policy is now testable without constructing solver or MPI runtime. This gives a clean seam to make future CLI migration deliberate. Next inspect whether the legacy ``d``/``r`` positional mode should become a named option, but first trace all invocations and preserve existing scripts.
- Suggested commit: ``refactor(cli): classify simulation invocation modes``.
### 2026-10-26 - Separate canonical CLI parsing from the legacy adapter

- Based on the latest inspected `hexin` commit `d5b96e999b35e3ae327affdbb62716d68f0889aa`; the user reports 723/723 tests passed on this baseline.
- `PrjCmdLineArgs.cpp` previously contained both the canonical simulation invocation/parser APIs and the compatibility conversion to `CmdLineOptions`. Moved the canonical implementations into `SimulationCommandLine.cpp`, leaving the old file with only the adapter.
- This keeps the generic parser (`codes/basic/include/CliParser.h`) independent, the simulation command-line API self-contained, and the legacy `Prj` surface at the compatibility edge. Registered the new source in each explicit test target that also compiles the adapter.
- Added a parity contract test proving the legacy adapter preserves debug mode, case order, repeated case paths, and the first-case convenience field. Moved parser, invocation-classification, and canonical-options tests into `simulation_command_line_test`, which compiles only the parser module and has no `Prj` or solver runtime dependency.
- Validation: direct MSVC compilation and execution passed 11/11 in the isolated `simulation_command_line_test` target and 41/41 in `prj_cmdline_test`. A fresh CMake configure completed and generated both Visual Studio targets; the generated projects list `SimulationCommandLine.cpp`. The subsequent CMake/MSBuild target build remained at `Checking Build System` without further output and was canceled, so CMake-driven compilation is unverified. The user reports 723/723 on the base revision; the full suite was not rerun.
- Architecture assessment: no runtime behavior or argv syntax changes. The test target still links the `Prj` compatibility adapter alongside the canonical parser, while parser-core tests remain independent of `Prj`.
- Suggested commit: `refactor(cli): separate parser from legacy adapter`.
### 2026-10-26 - Bring the Workbench GUI into the user and architecture manuals

- Based on current GitHub `hexin` HEAD `d5b96e999b35e3ae327affdbb62716d68f0889aa`; the user reports 723/723 tests on this baseline.
- Added the Chinese Sphinx user-guide entry `doc/source/workbench_user_guide_zh.rst`, added it to the main User Guides navigation, and linked it from `tools/oneflow_gui/README.md`. The module README remains the detailed source-adjacent reference; the Sphinx Workbench architecture page now separates implemented behavior from design gaps.
- Source review showed the Saved Workflows dialog does call `workflow_graph_runner.py`; the README statement that the GUI workflow editor was not connected was stale and has been corrected. The `.ofw` editor currently saves/loads/validates definitions but has no launch dispatcher. `ApplicationSettings.run_mode` is persisted and shown in Configuration but is not consumed by any launch path.
- Recorded both gaps in the user guide and architecture page. The immediate architecture choice is to wire the run-mode preference into an explicit dispatcher or remove the unused setting, then test that contract; later connect `.ofw` workspace modes to execution plans and run history.
- Validation: source inspection and documentation cross-link review. No GUI test suite or Sphinx build has been run in this slice; no full-suite result is claimed. The user's 723/723 result is for the starting revision.
- Suggested commit: `docs(workbench): add GUI user guide and current-state map`.

### 2026-10-26 - Connect the saved run mode to workflow draft selection

- Source inspection found that `WorkflowDialog` already accepts an `initial_mode` and creates a draft, but `WorkbenchHomeWindow` did not pass the persisted Configuration preference.
- Added a Qt-independent mapping from `multi_case` to `multi_directory` and from `dependency_workflow` to `case_graph`. Opening Saved Workflows now creates the matching draft; `single_case` keeps the saved-workflow list unchanged. This selects an editor draft only; it does not launch a run.
- Added three mapping tests and registered the application-settings test module in the GUI test runner. Stabilized an existing fallback test by using a guaranteed missing path rather than assuming the developer's default checkout path does not exist.
- Updated the GUI user guide, module README, architecture page, and AI handoff. The `.ofw` editor remains definition-only: integrating it with execution still needs a contract for artifact bindings and run-history records.
- Validation: Python syntax compilation succeeded; application-settings and workspace-persistence tests passed 17/17. The full C++/Python suite and Sphinx build were not run in this slice.
- Suggested commit: `refactor(gui): apply run mode to workflow draft selection`.

### 2026-10-26 - Run .ofw workspaces through the existing workflow runner

- Added `workspace_workflow_adapter.py`, a Qt-independent conversion boundary. Single, sequential, and independent workspace modes map to the existing one-process multi-directory workflow; dependency mode maps case IDs and edges to the existing graph runner.
- Added **Save and Run...** to the workspace editor. It saves the `.ofw` definition, opens a prefilled unsaved Saved Workflows draft, and lets the user confirm executable and paths before starting. Starting the run persists the workflow and run history through the existing repository.
- Dependency workspaces currently carry only ordering and failure-blocking edges. The `.ofw` format has no artifact input/output bindings, so it cannot hand generated files to consumers. The editor now states this limitation.
- Centralized `.ofw`-relative path resolution and serialization. Editing an existing file now preserves portable relative paths rather than resolving them against the process working directory.
- Added mode conversion and path round-trip tests. Validation: Python syntax compilation succeeded; the focused GUI/workflow suite passed 61 tests. The full repository suite and Sphinx build were not run.
- Suggested commit: `feat(gui): run saved workspaces through workflow runner`.

### 2026-10-26 - Add artifact bindings to .ofw dependency cases

- Bumped the `.ofw` schema to version 2 for optional `inputs` and `outputs`. Version-1 files without artifact fields migrate; version-1 files that contain those fields are rejected rather than silently reinterpreted.
- Added editable JSON binding columns to each case row. The workspace adapter forwards these values to the existing graph planner, which configures producer output settings, validates direct producer dependencies, stages outputs, and applies consumer settings.
- Reused `build_case_graph_plan()` as the canonical semantic validator for artifact IDs, input/output references, dependency edges, setting identifiers, and case-relative paths.
- Validation: Python syntax compilation succeeded; schema, conversion, persistence, GUI dialog, workflow plan, graph runner, staging, and process tests passed 83/83. The full repository suite and Sphinx build were not run.
- Suggested commit: `feat(gui): support artifact bindings in saved workspaces`.
### 2026-10-26 - Specify the ConfigLoader parser contract

- Based on `hexin` HEAD `de266a3b`; no runtime parameter or solver behavior was changed.
- Added `tests/io/ConfigLoader_test.cpp` and a focused test target covering scalar type/name/value parsing, declaration order, quoted values, multiline literal arrays, injected named-array length resolution, missing-resolver failure, and document replacement on repeated parse.
- Updated `doc/source/case_parameter_design.rst` to distinguish parser syntax guarantees from schema validation, runtime database application, and serialization guarantees.
- Validation: `ConfigLoader_test.exe` ran 5/5 tests, and `ctest --test-dir build -C Debug -R ConfigLoader --output-on-failure` passed all five after invoking the generated GTest discovery script. The CMake/MSBuild target returned failure with MSVC `C1902` PDB-manager mismatch errors; a clean target build is therefore unverified. A full Sphinx build of 109 sources passed using the installed Sphinx 9.2.0 with the built-in Alabaster theme override. The full C++ suite was not run.
- Remaining: resolve/verify the compiler toolchain path and ensure post-build GTest discovery runs automatically, then test type conversion at `ConfigDatabaseAdapter` and representative task resolution.
- Suggested commit: `test(config): specify legacy parameter parser contract`.


### 2026-10-30 - Isolate case I/O root setting selection

- Based on the latest inspected `hexin` revision before this slice, `eab1b83077e4ac5017e4c234f56250852971ecea`.
- Added `CaseIORootSettings` and `ReadCaseIORootSettings(const DataPara&)` in `SimuContextParameter.h`. `SetupCaseEnvironment()` now obtains the generic and domain-specific root settings through this adapter, then continues to apply them through the existing `Prj` APIs.
- Added focused tests for the complete root-setting snapshot, legacy `sourceGridCaseDir` fallback, and the requirement that an explicit `meshInputRoot` does not read an invalid legacy parameter.
- Updated the runtime architecture guide to distinguish parameter selection from global path application and file output.
- Validation: implementation and test changes are committed, but this environment did not build or run the target. The user's reported 963/963 result is the pre-change baseline and is not evidence for this slice.
- Architecture assessment: configuration lookup and legacy fallback selection now have a value-object boundary that can be tested without constructing the full solver runtime. This does not close R1: a real `SetupCaseEnvironment()` integration test still needs to exercise `ReadControlInfo(caseDir)`, `Prj` root application, `io_paths.txt`, and sequential-case reset together.


### 2026-10-30 - Pin project-root path resolution semantics

- Based on current `hexin` HEAD `4b063d1334298fc7802f71a5fcc12e989c173c22`; inspected `Prj::SetProjectBaseDir()`, `ResolveBaseDirectory()`, and existing `PrjIOPaths` tests.
- Added a focused contract test that makes the difference explicit: a relative `projectRoot` resolves against `Prj::current_dir`, while a generic output root set to `case` resolves against the active case directory.
- Documented the compatibility distinction in `doc/source/runtime_architecture.rst`. This is a behavior clarification, not a path-semantics change.
- Validation: source and test were edited through GitHub; no local build or CTest execution was available. The user-reported 967/967 pass remains the pre-change baseline and does not validate this test.
- System-level value: reduces ambiguity in the case configuration contract before extracting the production setup path. R1 remains partial until the real control reader, effective roots, diagnostics, sequential isolation, and failure cleanup are tested together.
