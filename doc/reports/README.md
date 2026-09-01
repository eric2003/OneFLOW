# Project reports

This directory contains maintained project-level reports, not per-run logs.

- Use stable filenames and update the existing report as the work evolves.
- Keep Markdown as the source of truth; regenerate the matching HTML after
  substantive changes.
- Store individual CI logs, scheduler records, and environment snapshots in
  GitHub Actions artifacts or the cluster's isolated run directory.
- Use `_drafts/` for local report drafts; that directory is ignored by Git.

Current maintained report:

- `gpu-backend-delivery.md`
- `gpu-backend-delivery.html`
- `one-dimensional-euler-hip-validation.md`
- `one-dimensional-euler-backend-optimization-plan.md`
- `one-dimensional-euler-performance-comparison-20260828.md`
- `one-dimensional-euler-performance-comparison-20260828.html`

The validation and optimization-plan reports are maintained as Markdown sources;
the dated performance comparison is maintained as both Markdown and standalone
HTML.
The validation report records the Kunshan correctness and performance evidence;
the optimization plan makes the stateful, persistent-device Euler backend the
current priority. WENO5, two-dimensional, and three-dimensional expansion are
deferred until the Euler acceptance gates pass.
The dated performance comparison report records the old trace and new stateful
benchmark methodology, normalized results, MPI resource comparisons, code changes,
build commands, and current claim boundaries.

The residual baseline database is maintained with the test definitions under
`test/baselines/`; it is not a per-run report.
