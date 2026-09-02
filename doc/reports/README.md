# OneFLOW project reports

This directory stores maintained project reports and selected historical snapshots.
Individual CI logs, Slurm records, and environment dumps belong in CI artifacts or
the isolated cluster run directories, not in this directory.

## Regression and residuals

- [Regression and residual integration (HTML)](regression/oneflow-residual-integration-20260901.html)
- [One-dimensional Euler validation](regression/oneflow-euler-validation.md)

## Performance

- [Current one-dimensional Euler performance report (Markdown)](performance/oneflow-euler-performance-current.md)
- [Current one-dimensional Euler performance report (HTML)](performance/oneflow-euler-performance-current.html)
- [2026-08-28 historical performance report (Markdown)](performance/oneflow-euler-performance-20260828.md)
- [2026-08-28 historical performance report (HTML)](performance/oneflow-euler-performance-20260828.html)

## Architecture and delivery

- [GPU backend delivery report (Markdown)](architecture/oneflow-gpu-backend-delivery.md)
- [GPU backend delivery report (HTML)](architecture/oneflow-gpu-backend-delivery.html)
- [One-dimensional Euler backend optimization plan](architecture/oneflow-euler-optimization-plan.md)

## Drafts

- _drafts/ contains local handoff notes and unfinished drafts. It is intentionally
  separate from maintained reports.

## Naming convention

Use oneflow-<topic>-<document-type>[-YYYYMMDD].<ext>:

- use a stable name such as current for the actively maintained report;
- keep the date for historical snapshots and point-in-time evidence;
- keep Markdown as the source of truth when a matching HTML report exists.

The residual baseline definitions remain under test/baselines/; they are test
inputs, not project reports.


## GitHub publication boundary

Suitable for the repository:

- maintained Markdown/HTML reports with conclusions, methodology, aggregate benchmark data, and reproducible code entry points;
- architecture notes, optimization plans, residual definitions, and sanitized resource tuples;
- historical snapshots when they explain the evolution of the implementation.

Keep outside the public repository:

- _drafts/, raw CI or Slurm logs, environment dumps, and temporary run directories;
- SSH keys, tokens, passwords, GitHub secrets, host-key material, or private cluster configuration;
- unredacted internal absolute paths, private hostnames, account identifiers, or scheduler/job metadata unless there is an explicit publication reason.

Public-facing cluster evidence should preserve the test meaning and resource tuple while replacing unnecessary internal identifiers with neutral labels. The _drafts/ directory is ignored by .gitignore and must remain a local working area.
