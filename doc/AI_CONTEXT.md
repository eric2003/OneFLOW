# AI and Developer Project Context

Last updated: 2026-11-01
Working branch: `hexin` only. Latest inspected baseline for this slice: `6e550260578a7ad3490793791f9a26feb1e670ae`. Current implementation HEAD before this handoff update: `7eff3d06099700cdd7026dca35cb843112b07332`.
User-reported local validation on the pre-change baseline: **1123/1123 tests passed**. This does not validate the shutdown-order change added afterward.

## Current follow-on: accelerator runtime shutdown ordering

The user reports **1123/1123 tests passed** on baseline hexin HEAD
6e550260578a7ad3490793791f9a26feb1e670ae. The current slice addresses a lifecycle
gap discovered while tracing SimuContext, EulerDomainStateRegistry, DeviceBuffer,
and AccelRuntime ownership.

- SimuContext owns EulerDomainStateRegistry; registered backend states may own
  device allocations.
- TeardownCase() already clears context states and releases the shared HIP flux
  backend before process finalization.
- FinalizeEnvironment() was independently callable but did not clear any
  remaining context-owned states before FinalizeAccelRuntime(). A caller that
  skipped TeardownCase() could therefore destroy backend state after runtime
  shutdown.
- FinalizeEnvironment() now clears the registry and releases HipFluxBackend's
  shared instance (when HIP is enabled) before FinalizeAccelRuntime(), then
  calls HXFinalize(). This is defensive process cleanup, not a replacement for
  full TeardownCase() of solver/global case state.
- A production-linked integration test inserts a state whose destructor
  records whether AccelRuntime is still initialized, then finalizes without
  TeardownCase() and checks destruction order.
- Design record: [Accelerator runtime shutdown order](plans/accelerator-runtime-shutdown-order.md).
- Runtime architecture guide updated in source/runtime_architecture.rst.

Industry references: C++ Core Guidelines R.1, C.31, C.36, C.37, E.4 and E.7.
The important distinction is that RAII manages each object's lifetime, while
the owning context must still order that lifetime before the allocator runtime
it depends on.

Status: **Implemented, unverified**. The user-reported 1123/1123 result predates
this slice. Build/run the production integration target, full Debug CTest,
and representative CPU regressions on the exact latest hexin HEAD. Also check
Windows/Linux CI. HIP-enabled compile/runtime coverage remains distinct from
the CPU ownership-ordering test.


## Previous mainline: R1 case configuration rollback and lifecycle

The mainline is to make the production case bootstrap exception-safe without confusing configuration rollback with full solver teardown.

Production path: `SimuImp::Run() -> RunCase() -> PreProcess() -> InitSimu() -> SimuContext::SetupCaseEnvironment() -> InitializeCaseConfiguration(caseDir) -> ReadControlInfo(caseDir) -> ReadCaseIOSettings() -> ApplyCaseIOSettings()`.

Current implementation:
- `codes/main/src/SimuContextCaseEnvironment.cpp` now wraps the production configuration bootstrap in `SimuCaseLifecycleDetail::ExecuteWithRollback()`.
- On bootstrap failure, `SimuContextDetail::ClearCaseConfigurationState()` clears partial DataField/DataPara state, the case log binding, and the project-root binding. The rollback helper preserves the original bootstrap exception if rollback itself throws.
- `IsEnvironmentReady()` remains false until the entire bootstrap succeeds.
- `SimuImp::RunCase()` still owns full `PostProcess() -> TeardownCase()` cleanup for solver/runtime resources on both setup and task failures. The new rollback is a narrower guarantee for callers of the public Context setup API; it does not replace full teardown.
- `tests/io/control_configuration_ingestion_test.cpp` now exercises the production Context setup path: the first script commits a value, the second script is missing, and the test verifies the partial value and project binding are cleared. It then creates the missing script and retries successfully on the same context.
- The separate two-case bootstrap test covers configured roots on one case and default-root rebinding on the next. It still uses the production configuration composition and real control reader.

System-level reason: case configuration now has a strong failure boundary independent of whether the caller is `SimuImp` or another Context consumer. This prevents partially ingested global configuration from leaking after a failed bootstrap. It does not yet prove full solver/device resource release through the actual `SimuImp::RunCase()` callback wiring.

## Validation status

- User reports 1096/1096 local tests passed on the pre-change baseline. No local build or CTest run is available in this environment for this slice.
- Runtime implementation commit: [6a325f7](https://github.com/eric2003/OneFLOW/commit/6a325f71ede48251df3d50e9e65d3a98f8b5dafe).
- Integration assertions and retry test: [4de3f86](https://github.com/eric2003/OneFLOW/commit/4de3f869affe58a5ac495a4cb90f12024bb3a6f7), [729811b](https://github.com/eric2003/OneFLOW/commit/729811b0996c929ec9b566bee643ab13b0897b75).
- Runtime architecture guide update: [abf2bd2](https://github.com/eric2003/OneFLOW/commit/abf2bd241d2139177630a2a99e7bb21d531e591a).
- Exact-HEAD Sphinx check succeeded; Windows and Linux CI were in progress at the last poll: [Sphinx](https://github.com/eric2003/OneFLOW/actions/runs/41280382266), [Windows](https://github.com/eric2003/OneFLOW/actions/runs/41280382252), [Linux](https://github.com/eric2003/OneFLOW/actions/runs/41280382250).
- This environment has no mounted OneFLOW Windows checkout, so no local VS2022/MSBuild/CTest result can be produced here.
- C++ source files remain ASCII-only by default; Python and mathematical presentation may retain useful Unicode characters.

## Remaining R1 gap

R1 remains **Partial**. The production integration target exercises the real `SimuImp::RunCase() -> PostProcess() -> TeardownCase()` wiring after a setup failure. A new success-path integration test is now implemented but unverified: it runs two separate `SimuImp` cases through the production `RunCasesWithProcessLifecycle()` orchestration, overrides `ToyModel` with a marker task, and checks case-state cleanup plus process-runtime continuity between cases. Remaining gaps are injected failures at individual solver/device release points, multi-rank MPI behavior, and representative full CFD regressions. The shared `OneFLOW_RuntimeObjects` CMake object collection now lets the executable and integration target consume the same production translation units without creating a static-library boundary. Do not substitute fake teardown and call it end-to-end.

The earlier investigation record remains `doc/plans/setup-case-environment-integration-test.md`. Its original findings are historical; this handoff and `doc/source/runtime_architecture.rst` record the newer rollback behavior.

Next single action: build and run the focused production lifecycle integration target on the exact latest commit, then run the full Debug CTest suite and representative grid/CFD regressions. If those pass, update the R1 evidence; keep teardown fault injection and multi-rank MPI as explicit remaining gaps. Detailed findings are in `doc/plans/setup-case-environment-integration-test.md`.
## Previous follow-on: accelerator buffer ownership

The latest branch focus has moved to accelerator resource reliability after the recent
transactional backend initialization work. Inspected starting HEAD: d9e25a7b671e6104c0d116ae5c561e5d12ba09d0 (Register accelerator runtime tests). The user reports 1119/1119 local tests passed on that pre-change baseline.

- DeviceBuffer<T>::Resize() previously leaked its replacement allocation if device-to-device copy failed.
- Its destructor and noexcept move assignment called a backend deallocator that may throw.
- The implementation now guards candidate allocations, preserves the old buffer on allocation/copy failure, rejects byte-count overflow and out-of-capacity host copies, commits replacement ownership before releasing the old allocation, and uses non-throwing cleanup from destruction/move assignment.
- Explicit Reset() reports release errors but clears object ownership first. Buffers must still be released before accelerator runtime finalization; a destructor cannot free through a backend that has already been finalized.
- Focused fake-backend tests cover failed growth cleanup/original-buffer preservation and non-throwing destruction.
- Design record: [DeviceBuffer resource lifecycle](plans/device-buffer-resource-lifecycle.md).

Implementation: [9728e46](https://github.com/eric2003/OneFLOW/commit/9728e46f97204a04a0bc188e770fde2e3d453e23). Focused test coverage and tracking backend: [df8dd00](https://github.com/eric2003/OneFLOW/commit/df8dd0008356e85ff0021cbb221535d9045069fd), [3bd2bf8](https://github.com/eric2003/OneFLOW/commit/3bd2bf84ca4f7b2abd71be94e2e77faf513ffffa), [e6d3e90](https://github.com/eric2003/OneFLOW/commit/e6d3e902e2f64dacb9ec1cb9a27cfbc3ca4cd2e5), [c19ef29](https://github.com/eric2003/OneFLOW/commit/c19ef29558f5a62ebe34e595d40976ad2af59aac), [a48c2b6](https://github.com/eric2003/OneFLOW/commit/a48c2b6f79651ec24e98d979a6a900d61a43dbf8), [74de922](https://github.com/eric2003/OneFLOW/commit/74de9222c90c7bce19602c209d4ddbdd1ab10b21).
Status at the time of that slice: implemented. The subsequent 1123/1123 local test report covers the DeviceBuffer changes, not the newer shutdown-order change. No local compiler/CTest execution is available in this environment. Next action is to build and run accel_runtime_test, then the full Debug CTest suite and CPU regression cases on the exact updated HEAD.

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

## 4. Current focus: R1 case configuration and lifecycle

Source and tests:
- `codes/main/src/SimuContextCaseEnvironment.cpp`
- `codes/main/src/SimuContextBootstrap.cpp`
- `codes/main/src/SimuContextCaseState.cpp`
- `codes/main/include/SimuCaseLifecycle.h`
- `tests/io/control_configuration_ingestion_test.cpp`
- `tests/main/simu_case_lifecycle_test.cpp`
- `tests/main/simu_imp_case_lifecycle_integration_test.cpp`
- `doc/source/runtime_architecture.rst`
- `doc/plans/runtime-refactoring-roadmap.md`

Verified contracts:
- Case setup publishes readiness only after the complete production bootstrap succeeds.
- A bootstrap failure rolls back partial DataField/DataPara state, case log binding, and project-root binding; the original setup exception is preserved.
- The same Context can retry after the failed manifest is repaired.
- `SimuImp::RunCase()` remains responsible for full case teardown; process finalization remains a separate owner.

Verified by production-linked integration: a real setup failure through `SimuImp::RunCase()` resets context-owned case state and leaves process runtime active until explicit finalization. Implemented, unverified: a successful two-case sequence using separate production `SimuImp` instances, checking per-case reset and process runtime continuity. Still open: injected teardown failures in the concrete release chain, multi-rank MPI, and full numerical regressions.
## 5. Next actions, in priority order

1. Check exact-HEAD Windows/Linux CI for the current `hexin` revision; distinguish queued, in-progress, and completed jobs.
2. Build and run the new successful sequential-case integration contract on the exact commit.
3. Run the full Debug CTest suite and representative grid/CFD regression cases.
4. If successful, update R1 evidence; then investigate teardown fault injection at concrete solver/device release points.
5. Keep multi-rank MPI validation as a separate acceptance item and reassess R1 status only against the required evidence.

Suggested next mainline: validate production sequential-case isolation, then close the remaining concrete cleanup-failure and MPI gaps—not another generic lifecycle helper or unrelated local refactor. The whole-project investment ranking is documented in `doc/source/architecture_priority_roadmap.rst`.
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


## 2026-11-01 continuation: shared production runtime and SimuImp integration

The user reports 1096/1096 local tests passed before this slice. This is a
pre-change baseline only.

Starting `hexin` HEAD: `0da2d68deeb7216c174e82630d75dab792103722`.

Implemented on `hexin`:
- `codes/CMakeLists.txt` now builds the recursively collected production
  source set as static `OneFLOW_Runtime`; `OneFlow.cpp` remains the sole
  executable entry point. The runtime target publishes the prior executable
  compile and link requirements so production consumers share one source graph.
- `SimuContext::IsProcessEnvironmentActive()` exposes a read-only process
  lifecycle observation.
- `simu_imp_case_lifecycle_integration_test` links the production runtime and
  invokes real `SimuImp::RunCase()` with a temporary case whose second control
  script is missing. It checks the full Context teardown through seeded task
  and solver-name state, configuration rollback, and process-runtime lifetime.
- Runtime architecture, lifecycle investigation, roadmap, and this handoff
  document the target boundary and acceptance gaps.

Important test contract: after the failed case, task resolution and expanded
solver names must be reset, partial configuration and project-root binding
must be absent, and process runtime must remain active until explicit
`FinalizeEnvironment()`.

No local build or CTest was run here. The 1096/1096 report is pre-change
baseline only. Exact-HEAD Windows/Linux CI and a fresh local build/test are
required because the executable was moved behind a new static-library target.

R1 is improved with real production `RunCase() -> PostProcess() ->
TeardownCase()` integration coverage, but remains Partial. Remaining:
individual teardown-step failure injection in the concrete chain, successful
repeated full SimuImp cases, multi-rank MPI, and CFD numerical regression.
Next action: inspect exact-HEAD CI and fix any build/link regressions before
expanding the integration suite.


Latest validation poll for `db7439f55c9b230004e09640f12f31fa860bcfad`:
- Windows CI queued: https://github.com/eric2003/OneFLOW/actions/runs/41281608533
- Linux CI queued: https://github.com/eric2003/OneFLOW/actions/runs/41281608515
- Earlier CI runs for the preceding C++ commits were still in progress; none had completed at this poll.
- No local CMake configure, build, or CTest was run in this environment. The user's 1096/1096 report remains the pre-change baseline.
