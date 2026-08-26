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

The one-dimensional Euler report is currently maintained as Markdown only. Its
cluster evidence is intentionally recorded without submitting new jobs during
the 2026-08-25 documentation pass; an HTML rendering can be added in a later
documentation-only update if needed.

The residual baseline database is maintained with the test definitions under
`test/baselines/`; it is not a per-run report.
