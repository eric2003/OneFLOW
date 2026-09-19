# AGENTS.md — guide for AI coding agents

This file is the entry point for AI coding agents working in this repository.
It indexes the maintained documents and states the rules that are not
negotiable. Read the linked documents before changing code or reporting
results.

## Documentation map

| Topic | Start here |
|---|---|
| Kunshan/DCU cluster workflow, verified environments, standard test suites | [`ci/kunshan/README.md`](ci/kunshan/README.md) |
| Report conventions, naming, publication boundary | [`doc/reports/README.md`](doc/reports/README.md) |
| Current 1D Euler CPU/DCU performance data | [`doc/reports/performance/oneflow-euler-performance-current.md`](doc/reports/performance/oneflow-euler-performance-current.md) |
| 2026-09-13 measurement method and baseline audit | [`doc/reports/performance/oneflow-euler-performance-20260913.md`](doc/reports/performance/oneflow-euler-performance-20260913.md) |
| CPU regression suite and residual baselines | [`test/README.md`](test/README.md), [`test/baselines/README.md`](test/baselines/README.md) |
| GoogleTest contract tests | [`tests/README.md`](tests/README.md) |
| Project state and handoff notes | [`doc/reports/architecture/oneflow-project-handoff-20260903.md`](doc/reports/architecture/oneflow-project-handoff-20260903.md) |

## Branch model

Contributions normally arrive from a fork, so a checkout can have two remotes:
`origin` (the fork) and `upstream` (this repository). Keep three kinds of
branches separate:

- **Baseline branch** (`master` in the usual layout) — kept identical to the
  upstream default branch and treated as a read-only PR baseline. Do not
  develop on it. Sync it with
  `git fetch upstream && git merge --ff-only upstream/master`.
- **Working branch** (`dev`, when the checkout has one) — the long-lived branch
  for daily work, in-progress features and notes. Pushing it to the fork is the
  backup step. Keep it current by **merging** the upstream default branch into
  it (`git merge upstream/master`); do not rebase it once it has been pushed,
  since that rewrites public history.
- **Topic branches** — short-lived, created from the upstream default branch
  only when a PR is explicitly authorized, by cherry-picking the function
  commits and excluding notes and unpublished docs. Delete them once the PR is
  done.

When working on a `dev` branch, read `doc/plans/oneflow-development-todo.md`
on that branch for the current working state — it lives on `dev` only.

## Working rules

1. **Regression before pull request.** Numerical-kernel or backend changes
   require the CPU five-case suite (normal `1e-8` and strict `1e-15`) and, for
   DCU/HIP changes, the HIP contract test. A change is not "validated" until
   the relevant suite passes and CI is green; do not open or report a PR as
   done on the basis of a build alone.
2. **Same-basis comparisons only.** `lifecycle_*_ms` in the benchmarks sums
   over `repeats`. Never divide datasets recorded with different `repeats`
   values, and state the basis whenever numbers are quoted.
3. **Markdown is the source of truth.** When a report has a matching HTML
   file, update both in the same change.
4. **Publication boundary.** Committed reports must not contain raw CI/Slurm
   logs, credentials, unredacted absolute paths, hostnames, account or job
   metadata. See `doc/reports/README.md`.
5. **Cluster work follows the standard suites.** Use the four suites defined
   in `ci/kunshan/README.md` (`cpu-regression`, `dcu-single`, `cpu-mpi`,
   `dcu-mpi`) with the documented workspace layout and resource tuples.
6. **Keep the port and the main solver separate.** The 1D Euler HIP code
   lives in `ports/kunshan/oneflow_1d_hip` and is built standalone with the
   DTK toolchain; the root build is CPU-only today. Do not claim a DCU
   capability from a root-build result.
