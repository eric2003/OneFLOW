# AI and Developer Project Context

Last updated: 2026-10-30
Working branch: `hexin`
Inspected starting HEAD for the current slice: `df567dba6f77abba725f05ce72cb32499046a653`.
User-reported pre-change baseline: 967/967 tests passed. This is not evidence for the changes below.

## Current focus: R1 case I/O root configuration contract

The active mainline is to make case setup configuration predictable and testable without coupling every contract test to the complete solver bootstrap.

Production path: `SimuImp::InitSimu() -> SimuContext::SetupCaseEnvironment() -> ReadControlInfo(caseDir) -> RequireGlobalDataBase().RequireDataPara() -> ReadCaseIOSettings() -> Prj root setters -> io_paths.txt`.

Current implementation:
- `codes/main/include/SimuContextParameter.h` defines `CaseIOSettings` and validates optional string entries before path application.
- `codes/main/src/SimuContextEnv.cpp` uses that snapshot and applies the resulting roots through the existing `Prj` compatibility APIs.
- An explicit `meshInputRoot` takes precedence over legacy `sourceGridCaseDir`; the legacy value is read only when mesh input is empty or `default`.
- `tests/main/simu_context_parameter_test.cpp` covers missing/valid/invalid optional strings, the complete root snapshot, and legacy fallback behavior. The current slice extends invalid-type coverage to every root key and `ioPathInfoFileName`; applicable settings are read before path-root mutation.

System-level reason: configuration lookup and type validation are now separated from mutation of global project paths. This reduces ambiguity in the runtime boundary, but does not prove the full setup lifecycle or per-case isolation.

## Validation status

- User-reported pre-change baseline: 967/967 tests passed.
- The updated type-contract test and documentation have not been built or run in this environment.
- No exact-commit CI result has been checked for this slice; do not claim the new test passes.
- Keep C++ source ASCII-only.

## Remaining R1 gap

A real integration test still needs to exercise the production `ReadControlInfo(caseDir)` path, effective `Prj` root application, `io_paths.txt`/database diagnostics, sequential-case reset, and cleanup after setup failure. The existing test targets intentionally avoid linking the full solver runtime; do not duplicate the parser or copy the executable's complete source list into a test target.

The existing investigation record is `doc/plans/setup-case-environment-integration-test.md`. It identifies the key constraints: the global database is statically initialized, but `ReadControlInfo` uses parallel coordination and database broadcast; `SimuContext::TeardownCase()` resets much more than path state; the production executable currently compiles a broad source list rather than exposing a reusable runtime library.

Next single action: determine whether a narrow production-owned case-setup orchestration boundary can preserve the real `ReadControlInfo(caseDir)` implementation while enabling setup success/failure and sequential-case tests. Do not mark R1 complete until effective path application and cleanup have test evidence.

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
- `codes/io/src/ParaFile.cpp`
- `tests/main/simu_context_parameter_test.cpp`
- `tests/main/prj_cmdline_test.cpp`
- `doc/plans/setup-case-environment-integration-test.md`

Verified behavior: missing optional root settings resolve to empty strings; present entries declared with a non-string type are rejected; generic and domain-specific roots are selected in a snapshot before being applied; explicit `meshInputRoot` wins over `sourceGridCaseDir`, and the legacy key is read only when mesh input is empty or `default`. Existing `Prj` tests cover path-root inheritance and reset rules.

Implemented, not yet verified for this slice: a test iterates over all eight root-setting keys and confirms that integer-valued entries produce a parameter-named type error. This strengthens the configuration schema contract, but is not a substitute for setup integration.

System-level concern improved: parameter interpretation is kept separate from mutation of process-global path state. Unchanged: real control-file loading through `ReadControlInfo(caseDir)`, effective `Prj` application, diagnostic file output, sequential-case isolation, and cleanup after setup failure.

## 5. Next actions, in priority order

1. Build/run `simu_context_parameter_test` and then the full suite on the current `hexin` revision; the reported 967/967 result is only the pre-change baseline.
2. Close the R1 integration gap with a narrow production-owned setup orchestration boundary. Preserve the real control reader; do not reimplement parsing or link a copied full executable source list into the unit test.
3. At the next milestone, refresh this handoff from the exact branch HEAD and inspect exact-revision CI evidence.

Suggested next mainline: a testable setup orchestration boundary, not another path resolver or a broad runtime rewrite.
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
