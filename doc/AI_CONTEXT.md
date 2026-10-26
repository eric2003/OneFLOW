# AI and Developer Project Context

Last updated: 2026-10-26
Working branch: `hexin`
Latest inspected local HEAD: `ea79bb66` (`origin/hexin` at inspection time).

The user reported 803/803 tests passing at the starting baseline. This working slice validates polygon-face connectivity before topology mutation and adds regression tests for undersized polygons and all-spans-before-staging behavior. The focused C++ target has **not** been rebuilt: MSVC reports C1902 PDB-manager errors across multiple geometry sources even after cleaning only the `OneFLOW_GeometryCore` target. A direct Clang override is incompatible with the generated project's `/scanDependencies` option. The geometry core outputs in this build tree therefore need a successful rebuild before treating the new tests as executed. Strict Sphinx HTML build passed for all 109 source documents with `-W --keep-going` after correcting two pre-existing unclosed emphasis markers on this page.

The goal of this slice is to make the input contract explicit: a polygon needs at least three node references, and the entire offset/connectivity span set must validate before any face lookup or staged topology mutation. Face identity is undirected/set-based, while the first-seen ring order remains the representative used by orientation and ownership processing. Repeated node IDs inside a ring remain outside this slice's validation contract.

This is a curated handoff snapshot, not a replacement for source code, tests, design pages, or CI. Recheck the branch HEAD and working tree before continuing.
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

## 4. Current focus: staged polygon-face validation

Current implementation files:

- `codes/geometry/src/FaceSolver.cpp`
- `tests/geometry/FaceSolverStageTest.cpp`
- `doc/source/face_stage_pipeline.rst`

The geometry path first validates the offset table size and every span; it now also rejects spans shorter than three references. Only after all spans pass does it construct rings, query the face lookup, and stage unique faces. The regression tests assert that malformed spans throw without leaving partial staged topology or desynchronized boundary-condition slots. The algorithm contract and rationale are recorded alongside the staged-face pipeline, including the distinction between undirected face identity and first-seen ring orientation.

System-level concern improved: malformed later input cannot leave earlier polygons partially staged. Unchanged: repeated node IDs within a polygon, MPI lifecycle, multi-case runtime, CLI parsing, and solver numerical kernels.

Validation evidence for this working slice:

- User-reported starting baseline: 803/803 tests passed; not rerun here.
- `git diff --check`: passed.
- C++ source ASCII check: passed for the changed implementation and test file.
- Strict Sphinx: passed for all 109 source documents with `-W --keep-going` after fixing two malformed emphasis markers.
- Focused C++ build/test: unverified due the C1902 errors above.
## 5. Next actions, in priority order

1. Restore a valid `OneFLOW_GeometryCore` build in the current build tree, then run the focused staging tests and full suite; do not cite the 803 baseline as verification of this change.
2. Inspect span arithmetic and malformed-connectivity behavior at the section reader boundary; keep any expanded validation contract paired with tests and docs.
3. At the next stable milestone, recheck the whole-system architecture risks in `source/architecture_review.rst` and update this handoff to the exact verified HEAD/results.

Suggested commit message: `fix(geometry): reject undersized polygon faces`.
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
