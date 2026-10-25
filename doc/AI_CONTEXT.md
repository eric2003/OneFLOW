# AI and Developer Project Context

Last updated: 2026-10-26
Working branch: `hexin`
Last inspected implementation commit: [`3c987064e4d547de4e10a07f261a3cebf2837407`](https://github.com/eric2003/OneFLOW/commit/3c987064e4d547de4e10a07f261a3cebf2837407)

This is the fast handoff page for an AI assistant or developer continuing work on OneFLOW. It is a curated status snapshot, not a replacement for source code, tests, design pages, or CI. Recheck the branch HEAD before editing: this page can lag behind the current repository revision.

## 1. Read this first

For a quick restart, follow this order:

1. **This page** — repository constraints, current focus, known gaps, and the next useful action.
2. [Runtime architecture and development guide](source/runtime_architecture.rst) — source-backed entry points, runtime lifecycle, configuration, and current boundaries.
3. [Runtime refactoring roadmap](plans/runtime-refactoring-roadmap.md) — staged plan, acceptance criteria, and per-slice findings.
4. [Documentation consolidation plan](plans/oneflow-documentation-consolidation-plan.md) — canonical documentation locations and maintenance rules.
5. The relevant focused tests and source files named below; then check CI for the exact commit being assessed.

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

## 3. Current focus: runtime case configuration and I/O roots

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

### Latest implementation slice recorded

Commit [`d8568fe91722af5d78d8083d6bae23d9446d952d`](https://github.com/eric2003/OneFLOW/commit/d8568fe91722af5d78d8083d6bae23d9446d952d) changes optional path-related string parameter reads in `SimuContextEnv.cpp`:

- Missing entries still return an empty string, preserving the existing fallback path.
- Present entries whose declared database type is not `HX_STRING` now raise `std::invalid_argument` rather than being read through an unchecked cast.
- The local helper is limited to environment setup; no public database API or generic path-resolver abstraction was added.
- The same rule is applied to `projectRoot`, `inputRoot`, `outputRoot`, `meshInputRoot`, `flowInputRoot`, `meshOutputRoot`, `flowOutputRoot`, `sourceGridCaseDir`, and `ioPathInfoFileName`.
- The follow-up documentation record is in [`3c987064e4d547de4e10a07f261a3cebf2837407`](https://github.com/eric2003/OneFLOW/commit/3c987064e4d547de4e10a07f261a3cebf2837407).

**This slice is not a completed R1 stage.** A focused automated test for the invalid-type setup path has not yet been added, and the prior record did not claim a local build or test run. Recheck the exact branch head and CI before relying on compilation or regression status.

## 4. Next actions, in priority order

1. Recheck current `hexin` HEAD and its Linux, Windows, and documentation CI runs.
2. Inspect the includes and type declarations used by `ReadOptionalStringParameter()`; add direct includes where required rather than relying on transitive headers.
3. Find a narrow, realistic test seam for the optional string parameter contract. Prefer a test of the actual helper or environment setup behavior over a test that only proves database type metadata.
4. Verify the path precedence and `sourceGridCaseDir` fallback contract with tests before extracting a resolver.
5. Update the roadmap and architecture guide with exact outcomes, test commands/results, CI links, and any remaining gap.
6. Only then decide whether R1 is complete and whether R2 (task selection and parameter semantics) should begin.

Do not mark an item complete because a code change was committed. Completion requires acceptance evidence.

## 5. Status vocabulary

Use these labels consistently in this page and related records:

- **Verified** — confirmed by source inspection, a named test, or a successful CI check; state which evidence.
- **Implemented, unverified** — code exists, but required test/build/CI evidence is not yet available.
- **Partial** — some behavior is covered, but an acceptance criterion remains open.
- **Planned** — not implemented.
- **Blocked** — progress depends on a named missing capability, test seam, or external condition.
- **Deferred** — intentionally postponed, with the reason recorded.

A green unrelated test suite is not evidence that a changed contract is correct. A queued or in-progress workflow is not a pass.

## 6. How to maintain this handoff

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
