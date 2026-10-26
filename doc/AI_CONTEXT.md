# AI and Developer Project Context

Last updated: 2026-10-26
Working branch: `hexin`
Latest inspected local HEAD: `299beb6f` (`origin/hexin` at inspection time).

The user reported 803/803 tests passing before commit `299beb6f`, which added a minimum-three-node check to polygon-face staging. This follow-up examines the topological meaning of that rule. OneFLOW's standard 2D fixed-node cell path derives two-node edges; the minimum-three rule is therefore specific to 3D polyhedral face loops. The CGNS adapter currently conflates NGON-only sections with polyhedral volume input, and does not import a 2D polygonal NGON section as cells. This is a documented adapter gap; do not generalize or remove the 3D face invariant to hide it.

This working slice narrows the staging API names to polyhedral faces, renames the dimension-aware ``IsBasicVolumeElementType`` predicate to ``IsBasicCellElementType``, adds a regression test for 2D triangle-to-edge cardinality, and documents cell/face dimension plus CGNS section roles. The full C++ target remains unverified: earlier MSVC attempts reported C1902 PDB-manager errors after cleaning only `OneFLOW_GeometryCore`, and the generated VS project rejected a Clang override because of `/scanDependencies`. The user's 803-test baseline is not verification of this slice. Strict Sphinx passed for all 109 source documents after these edits.

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

## 4. Current focus: dimension-aware mesh topology

Current implementation and design files:

- `codes/geometry/src/FaceSolver.cpp`
- `codes/geometry/src/GridElem.cpp`
- `codes/cgns/src/CgnsZsection.cpp`
- `codes/cgns/src/CgnsSection.cpp`
- `tests/geometry/FaceSolverStageTest.cpp`
- `doc/source/grid_topology.rst`
- `doc/source/face_stage_pipeline.rst`

Topology contract: a D-dimensional cell has codimension-one faces of dimension D-1. In 2D a linear edge face has two endpoint nodes; in 3D a polyhedral face is an ordered loop with at least three nodes. This is based on cell dimension, even for a surface mesh embedded in 3D coordinates.

Implemented slice: the generic-sounding polygon staging API is renamed to make its 3D polyhedral-face role explicit. Existing two-phase validation and minimum-three guard apply only to that role. A test records that TRI_3 creates three two-node faces. Current limitation: `HasPolygonSection` treats either `NGON_n` or `NFACE_n` as polyhedral routing evidence, while `SetElementTypeAndNode` accepts supported fixed-node cell types and does not convert 2D polygonal NGON cells into cells and edges. The classifier also reads global ``Dim::dimension`` rather than receiving cell dimension explicitly. Dimension-aware CGNS classification and 2D NGON import require a separate adapter design and acceptance tests.

System-level concern improved: cell dimension, face dimension, and connectivity role are now explicit in APIs/tests/docs. Unchanged: CGNS dispatch behavior, 2D NGON cell import, MPI lifecycle, CLI parsing, and solver numerical kernels.

Validation for this working slice: the user-reported baseline is 803/803 and was not rerun; `git diff --check` and strict Sphinx passed; C++ build/test is unverified due the known MSVC C1902 failure. C++ ASCII-only remains required.
## 5. Next actions, in priority order

1. Verify the renamed API and new 2D edge test with a working focused geometry build, then run the full test suite.
2. Design CGNS dispatch using base cell dimension and element-section roles; test NGON-only 2D polygon cells, 3D NGON plus NFACE polyhedra, and regular fixed-node meshes before changing routing.
3. At the next stable milestone, review the whole-system architecture risks and update this handoff with exact HEAD and validation evidence.

Suggested commit message: `refactor(geometry): clarify dimension-specific face staging`.
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
