# AI and Developer Project Context

Last updated: 2026-10-30
Working branch: `hexin`
Latest verified starting HEAD for this work slice: `4585de4e729e007e17e8bfbf6e2586ee7dc3c345`. Subsequent changes add production `ConfigDatabaseAdapter.cpp` coverage in `control_parameter_resolver_test`, and document its incremental ordered-commit semantics. The current source split places manifest orchestration in `ControlFileReader.cpp`, parsing and `ConfigDocument` commit in `ControlScriptParser.cpp`, database-backed array-size lookup in `ControlParameterResolver.cpp`, and broadcast/compression/decompression/diagnostic dumping in `ControlDatabaseIO.cpp`. `ParaFile.cpp` retains legacy parameter-inspection helpers; public declarations remain in `codes/io/include/ParaFile.h`. Recursive source discovery is expected to include the new translation units, but exact-revision build/CI evidence is still required.
User-reported pre-change local baseline for this work slice: 974/974 tests passed. This is not evidence for any later change.

## Current focus: R1 case I/O root configuration contract

The active mainline is to make case setup configuration predictable and testable without coupling every contract test to the complete solver bootstrap. The current intermediate test now pins the production ConfigDocument-to-DataBase conversion boundary, while the physical multi-file reader remains the next integration target.

Production path: `SimuImp::InitSimu() -> SimuContext::SetupCaseEnvironment() -> ReadControlInfo(caseDir) -> RequireGlobalDataBase().RequireDataPara() -> ReadCaseIOSettings() -> Prj root setters -> io_paths.txt`.

Current implementation:
- `codes/main/include/SimuContextParameter.h` defines `CaseIOSettings` and validates optional string entries before path application.
- `codes/main/src/SimuContextEnv.cpp` uses that snapshot and delegates root application and path diagnostics to `SimuContextDetail::ApplyCaseIOSettings()`, declared in `codes/main/include/SimuContextIO.h` and implemented in `codes/main/src/SimuContextIO.cpp`. The operation uses the existing `Prj` compatibility APIs. The implementation was moved out of the header so `Prj` and filesystem dependencies do not leak through an inline implementation.
- An explicit `meshInputRoot` takes precedence over legacy `sourceGridCaseDir`; the legacy value is read only when mesh input is empty or `default`.
- `tests/main/simu_context_parameter_test.cpp` covers missing/valid/invalid optional strings, the complete root snapshot, and legacy fallback behavior. The current slice extends invalid-type coverage to every root key and `ioPathInfoFileName`; applicable settings are read before path-root mutation.

System-level reason: configuration lookup and type validation are now separated from mutation of global project paths. This reduces ambiguity in the runtime boundary, but does not prove the full setup lifecycle or per-case isolation.

## Validation status

- User-reported pre-change baseline for the current work slice: 974/974 tests passed.
- The user-reported 974/974 result is the pre-change local baseline, not a result for the latest revision.
- Exact-HEAD CI for `d080aaa1b6a567ae3d0f2f5034f9332600c3e7ec`: [Windows](https://github.com/eric2003/OneFLOW/actions/runs/40954345525), [Linux](https://github.com/eric2003/OneFLOW/actions/runs/40954345454), and [documentation](https://github.com/eric2003/OneFLOW/actions/runs/40954345433) completed successfully. These results validate that revision only.
- The documentation workflow for `69a03cc27dd7ed7ab78c7516d4265138d00b20e8` completed successfully: [documentation check](https://github.com/eric2003/OneFLOW/actions/runs/41046094881). Windows and Linux workflows for that revision were still in progress at the latest check: [Windows](https://github.com/eric2003/OneFLOW/actions/runs/41046094753), [Linux](https://github.com/eric2003/OneFLOW/actions/runs/41046094882).
- No local test run is available in this environment. The user-reported 971/971 remains the pre-change baseline.
- The current Linux sandbox has CMake, GCC, Clang, Python, and Git, but no OneFLOW checkout is mounted and shell network/DNS access is unavailable. GitHub source edits are possible, but a complete local CMake/CTest run is not possible until a checkout/archive is available.
- Keep C++ source ASCII-only.

## Remaining R1 gap

A real integration test still needs to exercise the production `ReadControlInfo(caseDir)` path, effective `Prj` root application, `io_paths.txt`/database diagnostics, sequential-case reset, and cleanup after setup failure. Candidate fixtures confirmed in the current `hexin` tree are `test/grid/rae2822convert` (smaller multi-file control manifest) and `test/cfd/rae2822_roe_sa` (fuller CFD manifest; `simu.txt` declares `Solve`, which must not be executed by a setup-only test). The existing test targets intentionally avoid linking the full solver runtime; do not duplicate the parser or copy the executable's complete source list into a test target.

The existing investigation record is `doc/plans/setup-case-environment-integration-test.md`. It identifies the key constraints: the global database is statically initialized, but `ReadControlInfo` uses parallel coordination and database broadcast; `SimuContext::TeardownCase()` resets much more than path state; the production executable currently compiles a broad source list rather than exposing a reusable runtime library. The IO test target now covers the real document-to-DataBase adapter for typed scalar values and ordered duplicate replacement, but it still does not exercise the complete physical-manifest reader.

Cleanup ownership has now been distinguished explicitly in that investigation: `SimuImp::RunCase()` cleans up one case and leaves process runtime ownership to `SimuImp::Run()`, while the compatibility `SimuContext::SetupEnvironment()` rollback tears down both case and process runtime. Preserve these scopes; do not collapse them into a shared generic lifecycle. The post-`ReadControlInfo` application stage is now extracted as `SimuContextDetail::ApplyCaseIOSettings()` in `codes/main/include/SimuContextIO.h`. `SetupCaseEnvironment()` still calls the real reader first, validates a full snapshot, then invokes the extracted operation. The focused parameter test target now links the production `Prj` path sources and includes filesystem-backed diagnostic-path tests. This is source-level progress only; build and CTest evidence is still required.

Next single action: map the exact link dependency closure for `ReadControlInfo(caseDir)` and determine whether a small production configuration-runtime library can be shared by the executable and tests without duplicating CMake compile definitions or source lists. Then add a temporary multi-file integration fixture using the real reader, `SetUpParallelEnvironment()`, typed DataBase assertions, and `log/database.log` verification. Do not substitute a test parser or mock the database commit. Keep `tests/` for test code and `test/` for CFD cases. Do not mark R1 complete until effective path application, sequential-case isolation, and failure cleanup have test evidence.

## 1. Read this first

For a quick restart, follow this order:

1. **This page** — repository constraints, current focus, known gaps, and the next useful action.
2. [Whole-system architecture review](source/architecture_review.rst) — common CFD/CAE concerns, OneFLOW's cross-system risks, and the required zoom-out review at milestones.
3. [Runtime architecture and development guide](source/runtime_architecture.rst) — source-backed entry points, runtime lifecycle, configuration, and current boundaries.
4. [Runtime refactoring roadmap](plans/runtime-refactoring-roadmap.md) — staged plan, acceptance criteria, and per-slice findings.
5. [Documentation consolidation plan](plans/oneflow-documentation-consolidation-plan.md) — canonical documentation locations and maintenance rules.
6. The relevant focused tests and source files named below; then check CI for the exact commit being assessed.

Do not start by reading the whole repository. Use the architecture map to identify a concrete path, then inspect its callers, implementation, tests, and evidence.

## 2. Repository and collaboration constraints

- Repository: [`eric2003/OneFLOW`](https://github.com/eric2003/OneFLOW).
- Default working branch for this continuing effort: `hexin`. Base changes on the latest `hexin` HEAD, not `master`, unless the user explicitly changes this instruction.
- Keep changes merge-compatible with `hexin`; do not push or submit this work to `master`.
- Prefer small, evidence-backed changes that solve a real behavior or maintenance problem. Do not refactor merely to introduce a pattern.
- Follow direct call paths and state/resource ownership before changing interfaces.
- Preserve compatibility unless an intentional behavior change is stated and tested.
- Keep C++ source files ASCII-only.
- Do not claim a test ran, passed, or a CI job succeeded without evidence for the exact relevant revision.
- Keep public documentation aligned with current behavior. Distinguish implemented behavior, partial work, planned work, and unverified assumptions.
- The user values architecture discussion, meaningful tests, concise reasoning, and exact commit links.

## 3. Whole-system review checkpoint

The published [architecture review](source/architecture_review.rst) is now the common checklist for substantial AI-assisted changes. It deliberately separates general CFD/CAE concerns from verified OneFLOW behavior. Keep these cross-system risks visible while working locally: implicit global runtime state, configuration parsing versus validation/application, mutable mesh data and derived connectivity, workflow artifact provenance and execution locations, repeated/failed-run cleanup, and the gap between documented plans and actual test/CI evidence. At each meaningful milestone, state which system-level concern improved, what remains untouched, and whether the priority order should change.

## 4. Current focus: R1 case I/O root configuration

Source and tests:
- `codes/main/include/SimuContextParameter.h`
- `codes/main/src/SimuContextEnv.cpp`
- `codes/io/src/ControlFileReader.cpp`
- `codes/io/src/ControlDatabaseIO.cpp`
- `codes/io/src/ParaFile.cpp`
- `tests/main/simu_context_parameter_test.cpp`
- `tests/main/prj_cmdline_test.cpp`
- `doc/plans/setup-case-environment-integration-test.md`

Verified behavior: missing optional root settings resolve to empty strings; present entries declared with a non-string type are rejected; generic and domain-specific roots are selected in a snapshot before being applied; explicit `meshInputRoot` wins over `sourceGridCaseDir`, and the legacy key is read only when mesh input is empty or `default`. Existing `Prj` tests cover path-root inheritance and reset rules.

Implemented, not yet verified for this slice: a test iterates over all eight root-setting keys and confirms that integer-valued entries produce a parameter-named type error. This strengthens the configuration schema contract, but is not a substitute for setup integration.

System-level concern improved: parameter interpretation is kept separate from mutation of process-global path state. Unchanged: real control-file loading through `ReadControlInfo(caseDir)`, effective `Prj` application, diagnostic file output, sequential-case isolation, and cleanup after setup failure.

A new unverified test in `tests/main/prj_cmdline_test.cpp` pins the existing distinction that relative `projectRoot` resolves from `Prj::current_dir`, while ordinary relative I/O roots resolve from the case root. The runtime architecture guide records this compatibility behavior. This is source-backed but not yet built or run in this environment.

## 5. Next actions, in priority order

1. Check Windows and Linux CI for the exact source revision `a40941b661f53f54add332627cbfc939db057a34`; do not infer a pass from queued or in-progress runs.
2. Verify the recursive CMake source collector includes the parser and resolver translation units and inspect any CI link/build failures.
3. Preserve ConfigLoader's injected-resolver tests and the database-backed resolver tests as separate contracts.
4. Build a minimal target for the real `ReadControlInfo(caseDir)` path, initialized through `SetUpParallelEnvironment()`, with typed database assertions and `log/database.log` verification. Avoid copying the entire executable source list into a test target.
5. Extend from reader integration to effective `SetupCaseEnvironment()` path application, two-case isolation, and failure cleanup at the correct case/process lifecycle boundaries.
6. Reassess whether a reusable runtime library boundary is justified by this test and production dependency graph; do not introduce a generic lifecycle framework merely to make the test link.

Suggested next mainline: the real multi-file control-reader integration path, not another source split or path helper. The whole-project investment ranking is documented in `doc/source/architecture_priority_roadmap.rst`.
## 6. Status vocabulary

Use these labels consistently in this page and related records:

- **Verified** — confirmed by source inspection, a named test, or a successful CI check; state which evidence.
- **Implemented, unverified** — code exists, but required test/build/CI evidence is not yet available.
- **Partial** — some behavior is covered, but an acceptance criterion remains open.
- **Planned** — not implemented.
- **Blocked** — progress depends on a named missing capability, test seam, or external condition.
- **Deferred** — intentionally postponed, with the reason recorded.

A green unrelated test suite is not evidence that a changed contract is correct. A queued or in-progress workflow is not a pass.

## 7. Reusable stage records

For a substantial refactoring stage, use the [refactoring stage record template](plans/templates/refactoring-stage-record.md). Put detailed stage-specific problem statements, acceptance criteria, source paths, tests, decisions, and evidence in that record. Keep this page as a concise index of the overall current state and highest-priority next actions; link to detailed records rather than duplicating them.

## 8. How to maintain this handoff

Update this page when a milestone materially changes the state of the project, not after every command or minor edit.

For each meaningful update, record:

- inspected branch and base commit;
- the behavior or architecture question being addressed;
- current source path and relevant tests;
- what changed and what deliberately did not change;
- exact validation evidence, including failures and unavailable checks;
- unresolved gaps and the next recommended action;
- links to implementation and documentation commits.

Keep detailed reasoning and acceptance criteria in the relevant architecture page or roadmap. Keep reusable build/test commands in runbooks. Keep CI or runtime evidence in evidence records when a durable standalone record is useful. Avoid duplicating entire documents here; this page should point to them.

When starting a new AI session, first compare the current `hexin` HEAD with the implementation commit recorded above, inspect any changes since that commit, and then update this status snapshot only after verifying the new state.


## Evolutionary architecture

The new design page `doc/source/evolutionary_computation_architecture.rst` defines the long-term target: small computation units with explicit typed ports, preconditions/postconditions, side effects, lifecycle, failure policy, provenance, and verification evidence. The JSON Schema starter is `doc/specs/computation-unit.schema.json`. Both are proposals; no general plugin ABI or automatic pipeline runtime is claimed. The intended sequence is configuration ingestion first, then mesh/topology/boundary/result intermediate representations, and only later composable solver pipelines if real use cases justify them.


Latest source slice (implemented, unverified): `ReadOneFLOWScriptFile` overloads moved from `ParaFile.cpp` to `ControlScriptParser.cpp`; `GetParameterArraySize` moved to `ControlParameterResolver.cpp`. This makes the injected resolver seam in `ConfigLoader` explicit: parser tests can supply deterministic array-size resolution, while production uses the real database-backed resolver. The legacy `ParaFile.cpp` now retains inspection helpers. See `doc/plans/setup-case-environment-integration-test.md` for expected failure modes, partial-commit/transaction design question, and validation boundaries.

A focused `tests/io/control_parameter_resolver_test.cpp` and CMake target now cover literal array lengths, named lengths through the real DataBase, and missing-name failure. This is component evidence for the resolver only; full `ReadControlInfo(caseDir)` ingestion and parallel/database-dump integration remain open.


## Architecture investment and documentation policy

The ranked system-wide investment sequence is recorded in
`doc/source/architecture_priority_roadmap.rst`: (1) runtime lifecycle and
case isolation, (2) mesh core/topology/boundary lifecycle, (3) numerical
algorithm verification, (4) computation-unit contracts and dataflow
composition, (5) workflow/GUI, and (6) documentation automation. This is a
priority order, not a claim that lower-ranked work is unimportant.

Use source comments for stable API contracts and generate API reference pages
from Doxygen-compatible comments where practical. Comments should describe
intent, units, ownership/lifetime, invariants, preconditions, postconditions,
side effects, and failure semantics. Keep cross-file architecture, mathematical
derivations, numerical assumptions, experimental results, diagrams, and
reference comparisons in curated `doc/source` pages. Generated API docs
complement these guides; they do not replace them.

The architecture priority roadmap uses Sphinx math and text diagrams to make
the ranking, configuration transformation path, mesh representation layers,
and numerical verification chain explicit. Continue using equations, figures,
and diagrams where they explain actual algorithms or boundaries, and label
planned designs as proposals rather than implemented behavior.


## Documentation structure and source-comment policy

Keep `doc/plans/` separate from `doc/source/`. Plans record intended changes,
sequencing, acceptance criteria, risks, hypotheses, and unresolved questions.
The Sphinx source manual records current architecture, verified contracts,
algorithm derivations, scientific assumptions, figures, and stable design
knowledge. Promote completed plan conclusions into the relevant source page;
retain only useful decision history and remaining work in the plan. The policy
and directory map are described in `doc/README.md`.

Use Doxygen-compatible source comments for stable API contracts: purpose,
units, ownership/lifetime, invariants, preconditions, postconditions, side
effects, and failure semantics. Generate API reference pages from those
comments where practical. Keep cross-file architecture, derivations, formulas,
numerical verification, and conceptual diagrams in curated `doc/source`
pages. Source-generated API documentation complements, but does not replace,
the system and algorithm guides.

The 2026-10-30 IO adapter test adds two contracts: typed values reach the real
DataBase, and repeated declarations are committed in document order with the
later value replacing the earlier value. The adapter is incremental, not
transactional; changing that behavior requires separate compatibility analysis.
