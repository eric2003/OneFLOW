# OneFLOW Runtime Refactoring Roadmap

Status: In progress
Branch baseline: `hexin`
Last reviewed: 2026-10-26

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

**Status:** Selected for investigation; implementation and acceptance tests are not yet claimed complete.

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

**Status:** Planned; begin after R1 findings are recorded.

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

Append a short entry when a slice materially advances. Keep it factual and link the exact commit.

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
