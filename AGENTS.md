# AGENTS.md — guide for AI coding agents

This file is the entry point for AI coding agents working in this repository.
Before changing code, read [`doc/AI_CONTEXT.md`](doc/AI_CONTEXT.md) for the
current refactoring focus, verified findings, open gaps, and next actions.
Then follow the relevant architecture page, roadmap, tests, and evidence links.

## Fast documentation map

| Topic | Start here |
|---|---|
| Current AI handoff and ongoing refactoring state | [`doc/AI_CONTEXT.md`](doc/AI_CONTEXT.md) |
| Runtime entry points, lifecycle, configuration, and task dispatch | [`doc/source/runtime_architecture.rst`](doc/source/runtime_architecture.rst) |
| Runtime refactoring stages and acceptance criteria | [`doc/plans/runtime-refactoring-roadmap.md`](doc/plans/runtime-refactoring-roadmap.md) |
| Documentation ownership and consolidation | [`doc/plans/oneflow-documentation-consolidation-plan.md`](doc/plans/oneflow-documentation-consolidation-plan.md) |
| Kunshan/DCU cluster workflow and verified environments | [`ci/kunshan/README.md`](ci/kunshan/README.md) |
| Report conventions, naming, publication boundary | [`doc/reports/README.md`](doc/reports/README.md) |
| Current 1D Euler CPU/DCU performance data | [`doc/reports/performance/oneflow-euler-performance-current.md`](doc/reports/performance/oneflow-euler-performance-current.md) |
| 2026-09-13 measurement method and baseline audit | [`doc/reports/performance/oneflow-euler-performance-20260913.md`](doc/reports/performance/oneflow-euler-performance-20260913.md) |
| CPU regression suite and residual baselines | [`test/README.md`](test/README.md), [`test/baselines/README.md`](test/baselines/README.md) |
| GoogleTest contract tests | [`tests/README.md`](tests/README.md) |

## Branch and source-of-truth rules

- For the ongoing OneFLOW runtime, mesh, and Workbench refactoring, base work on
  the latest `hexin` HEAD unless the user explicitly gives a different target.
  Do not develop on, push to, or submit this work to `master`.
- Before editing, inspect the actual checked-out branch and current commit.
  The AI context page is a curated snapshot and may lag behind the repository.
- When records disagree, the current source and tests plus CI evidence for the
  exact commit determine implementation status. Architecture pages explain
  intended boundaries; plans and handoff notes are not proof that a feature
  or test exists.
- Keep changes compatible with `hexin` and avoid rewriting published history.

## Required working method

1. Trace the concrete entry point, direct callers, state/resource ownership,
   and failure/cleanup path before changing interfaces.
2. Prefer a focused change for a demonstrated behavior or maintenance problem.
   Do not refactor solely to introduce a pattern or a new abstraction.
3. Add or update tests for the changed contract at the narrowest meaningful
   boundary, followed by relevant integration or regression checks.
4. Report test commands, outcomes, failures, and platform limitations honestly.
   A build alone is not full validation; do not claim CI passed without checking
   the exact commit's workflow results.
5. Keep C++ source files ASCII-only. Preserve existing behavior unless an
   intentional behavior change is documented and tested.
6. At a meaningful milestone, update `doc/AI_CONTEXT.md` and the relevant
   architecture/roadmap entry with the verified base commit, change summary,
   evidence, remaining gaps, and next action. Do not create a command-by-command
   diary or duplicate full architecture documents.
7. Keep README and user/developer guidance consistent with current behavior.
   Clearly distinguish implemented, partial, planned, and unverified work.

## AI execution protocol for ongoing refactoring

Use this protocol for every continuation request. It is intended to be executable by
this assistant or another coding agent without relying on conversational memory.

### A. Establish the actual starting point

1. Inspect the current `hexin` HEAD, recent commits, working-tree/source state where
   available, and relevant CI results. Never assume that a handoff note or an old commit
   is still current.
2. Read `doc/AI_CONTEXT.md`, the relevant architecture/design page, the roadmap or stage
   record, and the directly related tests. Follow the real entry point and callers; do not
   scan or rewrite the whole repository without a reason.
3. Treat a user-reported test count as the pre-change baseline unless the user explicitly
   says otherwise. It is not evidence that new changes pass.

### B. Choose one high-value mainline

Before editing, identify one concrete mainline and explain briefly:
- the observed defect, risk, missing contract, or costly maintenance boundary;
- why addressing it now has greater system-level value than nearby alternatives;
- the exact behavior/acceptance criteria and the evidence that could prove or disprove it.

Prefer the change that closes a meaningful end-to-end gap, removes a real dependency or
ownership problem, or improves numerical confidence. Do not select work merely because it
is easy to rename or document. Keep the system-level goal visible while limiting the patch
to one coherent, independently reviewable slice. Do not start a second mainline in the same
turn unless the user explicitly asks.

### C. Make the smallest complete slice

1. State the intended behavior contract before changing implementation.
2. Inspect ownership, initialization, mutation, error propagation, and cleanup paths.
3. Add the smallest discriminating test at the narrowest useful layer; add integration,
   regression, or numerical checks when the contract crosses component boundaries.
4. Make the implementation change and update only the documentation that the change makes
   inaccurate or materially improves.
5. Preserve compatibility unless a behavior change is intentional, documented, and tested.
   Avoid broad rewrites and abstractions without a demonstrated need.

### D. Treat numerical knowledge as traceable engineering evidence

For a substantial numerical method or physical model, connect:
- governing problem, assumptions, units, symbols, sign conventions, and validity range;
- the formula actually implemented, including meaningful variants and limiting cases;
- original papers/books and suitable authoritative specifications or maintained manuals;
- OneFLOW source entry points and the difference, if any, between the reference and code;
- tests and cases that check invariants, analytic/manufactured solutions, conservation,
  convergence, benchmarks, or regression tolerances;
- known limits and evidence not yet collected.

Prefer primary papers, DOI/publisher pages, official project documentation, and stable
versioned sources. Other open-source implementations and public manuals are comparative
evidence, not automatic proof of correctness; record version/commit, relevant section or
file, and license. Use closed-source manuals only when lawfully available, cite edition and
section, and do not reproduce restricted code or prose. Do not add citations to trivial
helpers merely for appearance. A reference does not prove OneFLOW's implementation correct,
and a passing unit test does not prove the numerical model physically valid.

Use the loop: theory/reference -> OneFLOW formulation -> implementation -> focused tests
and representative cases -> discrepancy analysis -> updated explanation and evidence.
AI may help formulate hypotheses and propose tests, but must not silently change scientific
assumptions or treat its own generated output as independent validation. Changes to equations,
model variants, tolerances, and reference baselines require explicit review and reproducible
evidence.

### E. Validate and report honestly

- Run the most targeted feasible build/test first, then broader tests appropriate to the
  changed contract. Run representative CFD/grid cases when numerical or end-to-end behavior
  is affected.
- State exact commands and observed outcomes. Separate user-reported baseline, tests run
  locally, CI results for the exact commit, and checks that were not run or are blocked.
- Never infer a pass from compilation alone, an unrelated green suite, a queued workflow,
  or an old CI run.
- Use the status labels in `doc/AI_CONTEXT.md`: Verified, Implemented, unverified, Partial,
  Planned, Blocked, or Deferred, with evidence and remaining gaps.

### F. Keep commit history coherent

- Commit by logical work slice, not by file. Normally group implementation, its contract
  tests, and directly related architecture/roadmap/handoff updates into one commit.
- Use a second commit only when a change is independently useful, independently testable,
  intentionally reversible, or needs a distinct review boundary (for example, a separate
  fix discovered while validating the main slice).
- Do not create a separate commit for every file or every command. Do not combine unrelated
  changes just to reduce the commit count.
- Use concise messages describing the behavior or architectural outcome. Do not rewrite
  published history. Keep all work on `hexin` unless the user explicitly authorizes another
  branch.

### G. Leave a useful handoff

At a meaningful milestone, update `doc/AI_CONTEXT.md` and the relevant architecture/roadmap
record in the same logical change where practical. Record the inspected base/commit, chosen
mainline and rationale, files/contracts changed, evidence, remaining gaps, and the next
single recommended action. Keep the handoff concise; do not duplicate the full design guide
or write a command-by-command diary.

The goal is not to maximize code churn, test count, citations, or commit count. The goal is
to make OneFLOW easier to understand, test, validate, and improve from one evidence-backed
iteration to the next.

## Regression and publication rules

1. **Regression before pull request.** Numerical-kernel or backend changes
   require the CPU five-case suite (normal `1e-8` and strict `1e-15`) and,
   for DCU/HIP changes, the HIP contract test. A change is not validated until
   the relevant suite passes and CI is green.
2. **Same-basis comparisons only.** `lifecycle_*_ms` in benchmarks sums over
   `repeats`. Never divide datasets recorded with different `repeats` values,
   and state the basis whenever numbers are quoted.
3. **Markdown is the source of truth.** When a report has a matching HTML file,
   update both in the same change.
4. **Publication boundary.** Committed reports must not contain raw CI/Slurm
   logs, credentials, unredacted absolute paths, hostnames, account or job
   metadata. See `doc/reports/README.md`.
5. **Cluster work follows the standard suites.** Use the four suites defined in
   `ci/kunshan/README.md` (`cpu-regression`, `dcu-single`, `cpu-mpi`,
   `dcu-mpi`) with the documented workspace layout and resource tuples.
6. **Keep the port and the main solver separate.** The 1D Euler HIP code lives
   in `ports/kunshan/oneflow_1d_hip` and is built standalone with the DTK
   toolchain; the root build is CPU-only today. Do not claim a DCU capability
   from a root-build result.
7. **One fixed cluster workspace.** Kunshan work lives in a single
   `<workspace>/` root: `src/OneFLOW` is overwritten each round, `build/app` is
   the only incremental build, gate evidence goes to `evidence/<gate>/<jobid>/`,
   job logs go to `logs/`, and deprecated material goes to `_attic/<date>/`. Do
   not create per-job build/run directories or write logs to `$HOME`; see
   [`doc/runbooks/oneflow-kunshan-workspace-runbook.md`](doc/runbooks/oneflow-kunshan-workspace-runbook.md).
