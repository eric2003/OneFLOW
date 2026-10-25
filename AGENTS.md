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
