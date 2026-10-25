# AI and Developer Project Context

Last updated: 2026-10-26
Working branch: `hexin`
Latest R1 implementation slice: [`8cf20d87`](https://github.com/eric2003/OneFLOW/commit/8cf20d872d844b0fa68d99c5a59661be95ef0192) renames the legacy mesh-input fallback predicate to `ShouldUseLegacySourceGridCaseDir()` and adds tests proving that `case` and ordinary relative mesh roots do not select the legacy fallback. The user reported `100% tests passed out of 700` before this latest naming/test change; the new tests have not been independently built or run in this session. The runtime guide was aligned with checked optional-string lookup in [`da467a9c`](https://github.com/eric2003/OneFLOW/commit/da467a9c18cf5df0cc533410fd022d8dacc13afb), and its related-links line break was fixed in [`1cc6575d`](https://github.com/eric2003/OneFLOW/commit/1cc6575dc7baad91889b08a8339afc5230873802).

This is the fast handoff page for an AI assistant or developer continuing work on OneFLOW. It is a curated status snapshot, not a replacement for source code, tests, design pages, or CI. Recheck the branch HEAD before editing: this page can lag behind the current repository revision.

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

## 4. Current focus: runtime case configuration and I/O roots

The active refactoring stage is **R1: case configuration and I/O root resolution**. The goal is to clarify and test the existing contract before extracting a resolver or broadening the design.

### Source-backed execution path

`OneFlow.cpp -> Simulation::Run() -> SimuImp::InitSimu() -> SimuContext::SetupCaseEnvironment()`

Relevant implementation: `codes/main/src/SimuContextEnv.cpp`.

The setup path binds the active case, reads control parameters, applies project and generic I/O roots, applies mesh/flow-specific roots, handles the legacy `sourceGridCaseDir` mesh-input fallback, and writes path diagnostics.

### Established behavior from source inspection

- `Prj::SetPrjBaseDir()` resets project and domain-specific roots to prevent settings leaking between cases.
- Generic roots are applied before mesh/flow-specific roots.
- Empty or `default` domain roots inherit generic roots.
- Empty, `default`, or `case` generic roots resolve to the active case.
- The special `project` root requires a configured project root.
- Other relative roots resolve relative to the active case; absolute roots are used directly.
- `sourceGridCaseDir` is only a fallback when `meshInputRoot` is empty or `default`; an explicit mesh input root takes precedence.
- Missing `GetDataValue(name)` parameters throw. The underlying typed `DataObject` lookup uses a static cast, so it is not a safe runtime type check.
- Existing tests in `tests/main/prj_cmdline_test.cpp` cover many path/root rules. `tests/main/simu_context_test.cpp` primarily covers context state and task plumbing; the inspected tests did not directly exercise full `SetupCaseEnvironment()` precedence and legacy fallback.

See the runtime architecture guide and roadmap for the detailed evidence and limitations.

### Implementation and test seam

The current R1 implementation is split across these `hexin` commits:

- [`fd3889362408586e379d45d08a170abf2077efb9`](https://github.com/eric2003/OneFLOW/commit/fd3889362408586e379d45d08a170abf2077efb9) adds `SimuContextParameter.h`, a narrow lookup helper taking `const DataPara&`.
- [`6d72e31f93f55b00b622e122fa75eb689832694f`](https://github.com/eric2003/OneFLOW/commit/6d72e31f93f55b00b622e122fa75eb689832694f) routes production `SimuContextEnv.cpp` reads through that helper.
- [`1e438f84f600c689666f831ff0837ce3d5936e13`](https://github.com/eric2003/OneFLOW/commit/1e438f84f600c689666f831ff0837ce3d5936e13) adds focused tests for missing, valid string, and invalid-type parameters.
- [`a13a890bb31e2d84434fbe4297cf76aebc24c865`](https://github.com/eric2003/OneFLOW/commit/a13a890bb31e2d84434fbe4297cf76aebc24c865) registers the test target in CMake.

The helper preserves the optional-parameter contract:

- Missing entries still return an empty string, preserving the existing fallback path.
- Present entries whose declared database type is not `HX_STRING` now raise `std::invalid_argument` rather than being read through an unchecked cast.
- The local helper is limited to environment setup; no public database API or generic path-resolver abstraction was added.
- The same rule is applied to `projectRoot`, `inputRoot`, `outputRoot`, `meshInputRoot`, `flowInputRoot`, `meshOutputRoot`, `flowOutputRoot`, `sourceGridCaseDir`, and `ioPathInfoFileName`.
- The follow-up documentation record is in [`3c987064e4d547de4e10a07f261a3cebf2837407`](https://github.com/eric2003/OneFLOW/commit/3c987064e4d547de4e10a07f261a3cebf2837407).

**Status: Partial.** The user reported all 697 local tests passed on the code baseline before the latest fallback-policy change. The focused helper test covers missing, valid string, and invalid-type values; it now also tests whether empty/default mesh roots permit reading the legacy fallback and whether an explicit mesh root suppresses that fallback lookup. These latest additions still need the user's next local test run. The helper tests use a `DataPara` instance rather than the full `SetupCaseEnvironment()` lifecycle, so global-database binding and control-file loading remain integration-level acceptance criteria. The fallback predicate also has negative tests for `case` and an ordinary relative root. No GitHub CI pass has been confirmed for the latest revision.

## 5. Next actions, in priority order

1. Recheck current `hexin` HEAD and its Linux, Windows, and documentation CI runs.
2. Design the smallest integration test that exercises `SetupCaseEnvironment()` with a real temporary control file; first inspect `ReadControlInfo()` and global database initialization so the test does not accidentally duplicate runtime setup.
3. Assert explicit `meshInputRoot` overrides `sourceGridCaseDir`, and the legacy fallback is used only when mesh input is empty or `default`.
4. Verify the path reset between sequential cases and confirm output-path diagnostics remain compatible.
5. Re-run the full suite and record exact local/CI evidence after the integration test change.
6. Only after R1 is verified, decide whether R2 (task selection and parameter semantics) should begin.

Do not mark an item complete because a code change was committed. Completion requires acceptance evidence.

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
