# AI and Developer Project Context

Last updated: 2026-10-30
Working branch: `hexin`
Latest inspected starting HEAD for the 2026-10-30 continuation: `d00d86c64cf940f57106a68f62f406c8b6b36fc5`; latest source/documentation slice HEAD: `edacedd47f38c10b114f72c330af9243e16b5b8b`. Subsequent changes add production `ConfigDatabaseAdapter.cpp` coverage in `control_parameter_resolver_test`, and document its incremental ordered-commit semantics. The current source split places manifest orchestration in `ControlFileReader.cpp`, parsing and `ConfigDocument` commit in `ControlScriptParser.cpp`, database-backed array-size lookup in `ControlParameterResolver.cpp`, and broadcast/compression/decompression/diagnostic dumping in `ControlDatabaseIO.cpp`. `ParaFile.cpp` retains legacy parameter-inspection helpers; public declarations remain in `codes/io/include/ParaFile.h`. Recursive source discovery is expected to include the new translation units, but exact-revision build/CI evidence is still required.
User-reported pre-change local baseline for the 2026-10-30 continuation: 976/976 tests passed. This is not evidence for any later change.

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

## 2026-10-30 continuation after the 976/976 baseline

The latest inspected `hexin` starting HEAD was
`d00d86c64cf940f57106a68f62f406c8b6b36fc5`. The user reports 976/976 local
tests passed before this continuation; no new local test run is available here.
Windows and Linux CI for that exact HEAD were still in progress at inspection:
[Windows](https://github.com/eric2003/OneFLOW/actions/runs/41055632631),
[Linux](https://github.com/eric2003/OneFLOW/actions/runs/41055632649). Do not
report them as passing until their exact-run conclusions are checked.

The source review reconfirmed the production reader sequence:
manifest/script parsing -> process greeting/collection -> database broadcast
(one-process broadcast is a no-op after valid parallel setup) ->
`log/database.log` dump. A real serial integration test must call
`SetUpParallelEnvironment()`; static-zero parallel fields are not a valid
fixture. The smallest checked-in manifest fixture is
`test/grid/rae2822convert` and it selects the Grid task, which must not be
executed by a setup-only test.

CMake review found that `codes/CMakeLists.txt` recursively collects the
production code into the `OneFLOW` executable with global compile definitions
and dependencies, while IO tests select individual source files. No shared
configuration-runtime library currently exists. Do not copy the whole
executable source list into a test. Map the true dependency closure first and
create a shared production target only if a cohesive reusable boundary emerges.

Two documentation commits on `hexin` record this continuation:
- `cd806307af0566a889ba89c5da2baa3f1ae7d495` corrects the adapter source
  boundary: `ControlScriptParser.cpp` invokes the adapter implemented in
  `ConfigDatabaseAdapter.cpp`.
- `baea9f4e5342c84bde7e276b21a7f9ecefaeed45` records the reader execution
  sequence, serial initialization contract, CMake target boundary, and fixture.

Next action remains a real multi-file `ReadControlInfo(caseDir)` integration
test, following dependency closure mapping. R1 remains Partial; no evidence yet
for full reader integration, effective path application, sequential-case
isolation, or setup-failure cleanup.
## 2026-10-30 continuation: manifest discovery seam

Starting from the user-reported 976/976 local baseline, the latest `hexin`
source slice is `edacedd47f38c10b114f72c330af9243e16b5b8b`. No local build or
CTest run was performed in this environment. Exact-HEAD CI is queued:
[documentation](https://github.com/eric2003/OneFLOW/actions/runs/41057437914),
[Linux](https://github.com/eric2003/OneFLOW/actions/runs/41057437918),
[Windows](https://github.com/eric2003/OneFLOW/actions/runs/41057437955).

A narrow production extraction moved both `ReadScriptFileNameList` overloads
from `ControlFileReader.cpp` into `ControlManifestReader.cpp`. The public
contract in `ParaFile.h` now states that manifest paths are appended and their
order is preserved. `tests/io/control_manifest_reader_test.cpp` exercises the
case-directory overload using a temporary manifest and verifies blank-line
handling, case-root resolution, and order; it is registered with
`TextFileParser_test`. The new translation unit is also discovered by the
recursive production source collector.

This is component-level coverage for manifest discovery only. It does not
invoke script parsing/commit, `ReadMultiScriptFiles`, parallel coordination,
`ReadControlInfo(caseDir)`, or database dumping. R1 remains Partial. The
architecture rationale and exact boundary are in
`doc/plans/setup-case-environment-integration-test.md` and
`doc/source/case_parameter_design.rst`.

Next: inspect exact-HEAD Windows/Linux CI. If the new test target builds, map
the dependency closure from `ReadMultiScriptFiles` and
`ReadControlInfo(caseDir)`; then decide whether a production-shared
configuration-runtime target is justified. Do not broaden the test target by
copying the entire executable source list.


## 2026-10-30 continuation: manifest append contract

The user reports a new local baseline of 977/977 tests passed. This is the
pre-change baseline for the current continuation; it is not evidence for the
following GitHub edits.

Starting from `hexin` HEAD
`4cb14c0b28af7dd319cdc62e23075c8df1f066d9`, the manifest reader test was
strengthened to initialize its destination vector with an existing entry,
then assert that the entry remains at index zero while the three manifest
entries are appended in source order. This pins both parts of the public
contract: append semantics and stable manifest ordering. The production
reader is unchanged. Commit:
`13ce812c50b813b94afded7dbb1257a370876eaa`.

No local compiler or CTest run was available in this environment. The test
change must be validated on the user's Windows build and exact-HEAD CI before
being treated as verified.

System-level interpretation: manifest discovery is now a better specified
boundary, but it remains only the first stage of the R1 configuration path.
Do not expand this seam into a generic configuration framework. The next
high-value step remains mapping the real dependency closure for
`ReadMultiScriptFiles` and `ReadControlInfo(caseDir)`, then exercising
production parsing, typed database commit, serial coordination, diagnostic
output, sequential-case isolation, and failure cleanup.

Documentation direction: keep source comments focused on API contracts,
units, ownership, preconditions, postconditions, side effects, and stable
algorithm identifiers. Keep derivations, equations, assumptions, diagrams,
comparisons, and curated references in `doc/source`. Doxygen API output
should complement Sphinx rather than replace curated theory and architecture
pages. Prefer a small pilot on one stable module before adding Breathe or a
larger documentation toolchain.


## 2026-10-30 continuation: standardized engineering task specification

The user reports a pre-change local baseline of 977/977 tests passed. This is
user-reported evidence for the baseline, not a post-change result.

The task-expression requirements have been consolidated into maintained
specifications:
- `doc/AI_TASK_TEMPLATE.md`: canonical English engineering task specification.
- `doc/AI_TASK_TEMPLATE_zh_CN.md`: Chinese version for direct task preparation.
- `doc/README.md` and `doc/README_zh_CN.md`: documentation entry points.
- `doc/source/architecture_review.rst`: points to the task specification rather
  than maintaining a second, partially duplicated prompt.

The specification defines project context and exact base revision, observable
objectives, system-wide ROI prioritization, source-backed pre-change
investigation, implementation constraints, verification evidence levels,
documentation and scientific provenance, traceability, scope control,
acceptance criteria, and delivery records. The architecture review remains the
system-level priority framework; the task template governs individual task
execution; this file remains the current-state handoff.

R1 investigation reconfirmed that `ReadControlInfo(caseDir)` is a production
pipeline with ordered manifest discovery, multi-file parsing, typed database
commit, process coordination, database broadcast, and diagnostic output. A
valid serial fixture must initialize the existing parallel state explicitly;
static-zero fields are not a valid runtime contract. The existing focused
targets compile selected source files and do not yet provide a cohesive target
for the complete reader path. Do not copy the entire executable source list
into a test target. Continue mapping the dependency closure and identify the
smallest production-cohesive integration target.

Current implementation status remains:
- Manifest discovery has component-level coverage, including append and order
  semantics.
- R1 remains Partial: production multi-file reader integration, effective
  `SetupCaseEnvironment()` path application, sequential-case isolation, and
  failure cleanup remain unverified.
- No local build or CTest run was performed in this environment for the latest
  documentation and test changes. The 977/977 result is the user's pre-change
  baseline only.

Next action: inspect exact-HEAD CI for the test-contract change; then map the
production source/dependency closure for `ReadMultiScriptFiles` and
`ReadControlInfo(caseDir)`. Add an integration target only when its source
set and runtime initialization represent a real reusable boundary. Keep the
whole-system priority review active: runtime lifecycle first, mesh topology
and derived-data validity next, then numerical verification and documentation
automation as enabling work.


## 2026-10-30 continuation: ingestion dependency boundary and CI checkpoint

Starting revision: `4b2e5d496ed359c1e3c3645413c4541dc81122f7` on `hexin`. The user reports 977/977 local tests passed before this continuation; this is a pre-change baseline only.

The source-level dependency map and revised test layering are recorded in [the R1 integration investigation](plans/setup-case-environment-integration-test.md#2026-10-30-follow-up-source-level-dependency-graph-checkpoint), committed as `10a21b46b3f071f233dcae854d24dad09c1c91c2`. The key decision is to test the real configuration-ingestion composition first (manifest -> multiple production script parsers -> typed DataBase commit), separately from the broader `ReadControlInfo(caseDir)` orchestration (parallel coordination, broadcast, logging, PIO and database dump). This avoids turning the first test into an improvised miniature of the full executable.

Exact-HEAD CI observed for starting revision `4b2e5d496ed359c1e3c3645413c4541dc81122f7`:
- Sphinx build: **success** — [job](https://github.com/eric2003/OneFLOW/actions/runs/41059549136/job/124075558246).
- Linux OneFLOW Actions: **success** — [job](https://github.com/eric2003/OneFLOW/actions/runs/41059549091/job/124075553594).
- Windows OneFLOW Actions: **in progress** at the time of inspection — [job](https://github.com/eric2003/OneFLOW/actions/runs/41059549210/job/124075554479).

No source code or tests were changed in this continuation; no local build or CTest was run. The documentation-only dependency-map update is now at `10a21b46b3f071f233dcae854d24dad09c1c91c2`. R1 remains **Partial**. Next: implement the real multi-file configuration-ingestion integration target with RAII-managed temporary fixtures, actual DataBase commits, typed assertions, and explicit failure semantics. Then add the broader `ReadControlInfo(caseDir)` orchestration test with a correctly initialized serial/MPI runtime and diagnostic output. Recheck CI on the new exact HEAD before claiming validation.


## 2026-10-30 continuation: production configuration-ingestion composition

Starting `hexin` revision: `4b359903f4ffbf3bf6c595f066d920bf6b26e5b3`. User-reported local baseline: 977/977 tests passed before this change; do not present it as post-change validation.

Implemented on `hexin`:
- Moved `ReadMultiScriptFiles` from `ControlFileReader.cpp` into `ControlScriptParser.cpp`, separating ordered parsing/commit from runtime parallel orchestration and database dumping.
- Changed its input to `const std::vector<std::string>&` and documented ordering plus per-file incremental commit semantics in `ParaFile.h`.
- Added `tests/io/control_configuration_ingestion_test.cpp`: a unique temporary case with a real manifest and two scripts is read through the production manifest reader and multi-file parser, then typed values and later-file replacement are asserted against the actual global DataBase.
- Expanded and renamed the combined CMake target to `control_configuration_test`; existing resolver/adapter tests remain included. The target links the production parser, manifest, resolver, adapter, database and low-level support files, but deliberately excludes `ControlFileReader.cpp` and its process/broadcast/dump dependencies.
- Detailed implementation scope and limitations are recorded in `doc/plans/setup-case-environment-integration-test.md` (commit `8ef8c740fbb27e85d1ea5eb9f83950da8a2097b6`).

Validation status: **Implemented, unverified**. No local compiler/CTest run was available here. The 977/977 result is the pre-change baseline. This test does not cover `ReadControlInfo(caseDir)`, parallel initialization, broadcast, database diagnostic output, `SetupCaseEnvironment`, sequential-case isolation, or failure cleanup. R1 remains **Partial**. Next: inspect exact-HEAD Windows/Linux CI, resolve any CMake/link errors, request and record the user's post-change local test result, then add a separate faithful orchestration test only after its runtime dependency closure is explicit.


## 2026-10-30 continuation: CI regression diagnosis and CMake correction

Starting revision: `326a576e9271ed00efaca30a769d70b87f18ebe5`. User-reported local result remains 977/977 from before the changes; it is not post-change evidence.

Exact-HEAD CI for `326a576e...` exposed the same configuration failure on both platforms, before compilation or tests:
- Windows configure failed at `tests/io/CMakeLists.txt:59`: CMake could not find `control_configuration_test.cpp`, then reported no sources for `control_configuration_test`.
- Linux configure failed for the identical reason.
- Sphinx build passed.

The source inventory confirms the existing resolver/adapter test file is still named `tests/io/control_parameter_resolver_test.cpp`; the target had been renamed but its source path had not. Corrected the CMake source entry to the existing filename in commit `d9180818eed6946130cfdd29a96477ff8c932e7d`. No production C++ behavior was changed by this correction.

At the time this handoff was written, Windows and Linux CI for `d9180818...` were both in progress:
- [Windows job](https://github.com/eric2003/OneFLOW/actions/runs/41063294279/job/124088158630)
- [Linux job](https://github.com/eric2003/OneFLOW/actions/runs/41063294280/job/124088158905)

The new integration test currently verifies the successful two-file manifest -> production parser -> typed global DataBase path. It does not yet verify missing-file failure semantics, `ReadControlInfo(caseDir)`, MPI/serial orchestration, diagnostic output, case path application, sequential-case isolation, or cleanup at the real case runtime boundary. R1 remains **Partial**. Next action: wait for exact-HEAD CI, inspect any subsequent compile/link/test failures, and only then extend the integration behavior. Avoid treating a successful CMake configure as proof that the expanded source target compiles or runs.


## 2026-10-30 continuation: const DataObject lookup compile fix

User-reported MSVC diagnostic while compiling `control_configuration_ingestion_test.cpp`: `DataBase.h:152` could not resolve `GetDataValue<T>(&data)` when `data` was `const DataObject&`.

Root cause verified in the declarations:
- `DataObject*` overload declared `int iElement = 0`.
- `const DataObject*` overload declared `int iElement` with no default, although its definition takes the same two parameters.
- The const-DataBase `GetDataValue` implementation passes only the `const DataObject*`, so overload resolution fails.

Corrected `codes/database/include/DataObject.h` by giving the const overload the same default element index (`= 0`). This preserves const correctness and aligns the two scalar lookup APIs; it does not cast away constness or change runtime behavior. Commit: `e8fa80ec66ed78790e7e9b13b62e22503823cf94`.

The ingestion integration test's `const DataBase&` lookup now exercises this const overload directly. This is useful API-contract coverage: const database access must support scalar element-zero reads just as mutable database access does.

Validation: source-level fix only at this point. No local MSVC build or CTest has been run by the assistant. The 977/977 result remains the pre-change baseline. Next: verify exact-HEAD Windows/Linux CI; inspect any compile/link/runtime failures before expanding R1. R1 remains **Partial**.


## 2026-10-30 continuation: database contracts and documentation

User-reported local validation before this continuation: 978/978 passed. The previous exact HEAD `cc73e947ceb0c9b9b88477f81aa0578012ba827b` had successful Windows, Linux, and Sphinx checks.

Changes made:
- Added API-level documentation to `DataObject.h` for the type-erased value boundary, borrowed raw storage access, and the preconditions of mutable/const `GetDataValue<T>` (non-null object, exact stored type, valid index). The two overloads have matching default index zero.
- Documented `DataEntry`/`DataPara` ownership, borrowed lookup lifetimes, and same-key update behavior in `DataPara.h`.
- Documented `DataBase` ownership of parameter and field registries, and mutable/const named lookup contracts in `DataBase.h`.
- Extended the real configuration-ingestion test to exercise both mutable and const `DataBase` scalar reads.
- Added `doc/source/database_architecture.rst` and linked it from the Sphinx index. It explains ownership, type erasure, parameter replacement, constness, field pointer registration, serialization, and the boundary between ingestion integration and full runtime integration.

Commits in this documentation/API slice: `e1e3a231226635709156b939f686f47661268ffb`, `8504f6a5c2e85efe146823216f6b19d4ed7fef03`, `afe134cc5a5451a94341fb7b05f458aebbcb3999`, `c4b816cbd865223360f685a23ccd96dd4344d993`, `a6abebde9dc3120240ef2edf713b14a2c8ea8bc5`, `eb8590aef26e2ac0a47c9cb4ed85aa5ef5ca48ed`.

Validation status: these changes have not yet been validated by a new local run or exact-HEAD CI. The 978/978 result is the pre-change baseline for this slice. R1 remains Partial. Next inspect exact-HEAD CI, then continue with explicit missing-file/parse-failure semantics and the separate runtime orchestration layer. Avoid implying the database guide proves unsupported behavior; it records current source-backed contracts and clearly separates remaining gaps.


## 2026-10-30 continuation: configuration-ingestion failure semantics

The user reports a pre-change local baseline of 978/978 tests passed. Do not
treat this as validation of the following commits.

Starting `hexin` HEAD was
`3fadcbacafee01081b65bc889218045566dc6a46`. Sphinx passed on that exact
revision; Windows and Linux build/test jobs were still in progress when
inspected.

Source review of `ConfigLoader`, `ControlScriptParser`,
`ConfigDatabaseAdapter`, and `Prj::OpenFile` established the following
current behavior:

- The multi-file reader is sequential and incremental across files.
- Each file is fully parsed into a fresh `ConfigDocument` before its adapter
  commit begins.
- Missing-file failure occurs during file open and leaves prior file commits
  intact.
- A parse-time failure prevents the current file's document from reaching the
  adapter, so declarations parsed earlier in that same file are not committed.
- The adapter's commit is itself not atomic across entries; a later entry
  failure can leave earlier entries from the same document committed.

Added tests for the first two failure boundaries in
`tests/io/control_configuration_ingestion_test.cpp`, committed as
`7c3d488665a35cd5aa293822c64ce9f8b5b91cd3`. The investigation plan was
updated in `doc/plans/setup-case-environment-integration-test.md`, committed
as `df4391a13de37950b24e3889df36eeb7dec5b65b`.

The tests intentionally specify current behavior; they do not introduce
whole-manifest rollback. No local build or CTest was run in this environment.
The user's 978/978 result is the pre-change baseline only. Recheck exact-HEAD
Windows, Linux, and Sphinx CI after the latest documentation commit.

R1 remains Partial. The next larger boundary is still the real
`ReadControlInfo(caseDir)` orchestration, including initialized serial
parallel state, broadcast path, and `log/database.log`; after that, verify
`SetupCaseEnvironment` path application, sequential-case isolation, and
failure cleanup. Do not expand the test target by copying the full executable
source list.


## 2026-10-30 latest continuation: script reader boundary

User-reported local baseline before this slice: 981/981 passed. The starting
`hexin` HEAD was `0e0054d658db092e1cfbc54195ac87d6ef3caf84`; do not treat
that baseline as validation of this slice.

The case-level ingestion test now calls the production
`ReadPrjScript(caseDir)` composition rather than manually invoking manifest
discovery and multi-file parsing. To keep this production entry point
independently linkable, `ReadPrjScript` overloads were moved from
`ControlFileReader.cpp` into the new
`codes/io/src/ControlScriptReader.cpp`. `ControlFileReader.cpp` now retains
the `ReadControlInfo` orchestration boundary (parallel greeting, broadcast,
and database dump). The new source is included by recursive production source
collection and explicitly by `control_configuration_test`.

This separation is an intentional dependency boundary: script ingestion is
distinct from process orchestration. It does not claim full
`ReadControlInfo(caseDir)` coverage. Exact-HEAD Windows/Linux CI and local
CTest are required; no local build was run here. The 981/981 report predates
these edits.

Next: inspect the exact link dependency closure of `ReadControlInfo` and
consider a cohesive production target only if it is useful to both executable
and tests. Preserve real serial/MPI initialization, DataBook broadcast,
logging, PIO, and database dump semantics. R1 remains Partial until
`ReadControlInfo`, effective Case path application, sequential-case
isolation, and failure cleanup are verified.


## 2026-10-30 continuation: production reader orchestration slice

The user reports 981/981 tests passed before this slice. This is a pre-change
baseline only. Starting `hexin` revision:
`144cde3f472fc73d2a8633f84b5eaea874308de8`.

Implemented source changes:
- Corrected the case database dump to use `PIO::OpenCaseFile` rather than
  calling `Prj::OpenCaseFile` directly. The latter performs path resolution
  and file creation but no parallel-rank filtering; the old case dump could
  therefore open/truncate the same output on every MPI rank while its close
  operation was owner-filtered.
- Split explicit path-based PIO operations into
  `codes/io/src/PIOFileAccess.cpp`. The remaining `PIO.cpp` operations are
  task-state-bound convenience entry points. This is a real dependency
  boundary, not a duplicate test implementation.
- Extended `control_configuration_test` to compile the production
  `ControlFileReader.cpp`, `ControlDatabaseIO.cpp`, `PIOFileAccess.cpp`,
  and `LogFile.cpp`.
- Added a serial integration test that initializes the production serial
  parallel state, invokes `ReadControlInfo(caseDir)`, and checks manifest
  ordering/override, typed global DataBase values, `log/database.log`
  contents, and the case-bound diagnostic log.

The new test covers the real serial reader orchestration and the serial
broadcast short circuit. It does not establish MPI collective correctness or
multi-rank file ownership at runtime. Full `SetupCaseEnvironment()`,
effective path application through that entry point, sequential-case
isolation, and exception cleanup remain unverified. R1 remains Partial.

The source-backed dependency and ownership analysis is recorded in
`doc/plans/setup-case-environment-integration-test.md`.

Validation: no local compiler or CTest run was performed here. The user's
981/981 result predates this slice. Inspect exact-HEAD Windows/Linux CI and
request a fresh local test result before claiming the new integration test
passes.


## 2026-10-30 continuation: geometry-core PIO link dependency

The user reported MSVC linker errors from `GridSummaryAdapterTest.exe`:
`OneFLOW_GeometryCore.lib(Zone.obj)` and `GridGroup.obj` referenced
explicit-path PIO methods that were no longer defined in `PIO.cpp`.

Root cause: `tests/geometry/CMakeLists.txt` builds the shared
`OneFLOW_GeometryCore` static library from an explicit source list. After
explicit-path PIO methods were moved to `codes/io/src/PIOFileAccess.cpp`,
this target still listed only `PIO.cpp`. The production executable's
recursive source collection was not affected, and the IO ingestion test
already listed the new source; the geometry test core was the missing
consumer.

Added `PIOFileAccess.cpp` beside `PIO.cpp` in the geometry-core source list.
This fixes the shared library dependency for all geometry tests linking
`OneFLOW_GeometryCore`, rather than patching only `GridSummaryAdapterTest`.

Commit: `eab78cbcc4a0b9ee1fb59a439844a67c6bf7a687`.
No local build or CTest was run here. Reconfigure/rebuild the Windows Debug
geometry test target and run the full suite; this commit is not yet validated
by the user's local test result. Also inspect exact-HEAD CI before treating
the branch as green.


## 2026-10-30 continuation: Case lifecycle integration boundary

The user reports 982/982 local tests passed before this continuation. A fresh
source review mapped the real production path across
`Simulation::RunImpl`, `SimuImp::RunCase`, `SimuContext::SetupCaseEnvironment`,
`TeardownCase`, and process finalization.

Important ownership result: one process context owns parallel/accelerator
initialization and finalization; each case gets a distinct `SimuImp` /
`SimuContext`; case setup and teardown run inside the case cleanup boundary.
A failure in `SetupCaseEnvironment` still reaches `PostProcess` /
`TeardownCase`, after which the outer process lifecycle finalizes the shared
runtime.

The IO path reset contract already has a component test:
`PrjIOPaths.DomainRootsResetWhenBindingNextCase`. It configures project and
domain-specific roots on one case, rebinds a second case, and verifies all
mesh/flow roots resolve under the second case. This proves the `Prj` path
component reset, not the complete `SetupCaseEnvironment` composition.

The remaining R1 gap is therefore not missing isolated tests for each
operation: the real `ReadControlInfo`, IO settings application, and
`Prj` root reset are individually exercised. What remains untested is their
production composition inside `SimuContext::SetupCaseEnvironment`, plus
cleanup after failure at that boundary. The investigation and next extraction
proposal are recorded in
`doc/plans/setup-case-environment-integration-test.md`.

No source code or post-982 test result is claimed in this slice. The latest
hexin revision for this note is recorded by its commit history; check exact
HEAD CI before declaring the preceding PIO source-list fix green.
