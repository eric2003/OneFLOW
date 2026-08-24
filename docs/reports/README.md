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
