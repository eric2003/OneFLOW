# OneFLOW project reports

This directory stores maintained project reports and selected historical snapshots.
Plans, handoffs, runbooks, and validation evidence live in their dedicated sibling
directories under `doc/`; they do not belong here. Individual CI logs, Slurm records,
and environment dumps belong in CI artifacts or isolated cluster run directories.

## Regression and residuals

- [Regression and residual integration (HTML)](regression/oneflow-residual-integration-20260901.html)
- [One-dimensional Euler validation](regression/oneflow-euler-validation.md)
- [One-dimensional Euler validation (HTML)](regression/oneflow-euler-validation.html)

## Performance

- [2026-10-14 standalone 1D Euler MPI/DCU strong-scaling benchmark](performance/oneflow-1d-euler-mpi-dcu-benchmark-20261014.md)
- [2026-10-19 standalone 1D Euler rank/DCU matrix](performance/oneflow-1d-euler-mpi-rank-dcu-matrix-20261019.md)
- [2026-09-13 CPU/DCU performance comparison](performance/oneflow-euler-performance-20260913.md)
- [Current one-dimensional Euler performance report](performance/oneflow-euler-performance-current.md)
- [Current one-dimensional Euler performance report (HTML)](performance/oneflow-euler-performance-current.html)
- [2026-08-28 historical performance report](performance/oneflow-euler-performance-20260828.md)
- [2026-08-28 historical performance report (HTML)](performance/oneflow-euler-performance-20260828.html)

## Delivery

- [GPU backend delivery report](delivery/oneflow-gpu-backend-delivery.md)
- [GPU backend delivery report (HTML)](delivery/oneflow-gpu-backend-delivery.html)

## Related project documents

- [Architecture](../architecture/)
- [Plans and living TODO](../plans/)
- [Handoff](../handoff/)
- [Runbooks](../runbooks/)
- [Validation evidence](../evidence/)
- [Drafts](../_drafts/)

## Naming convention

The authoritative repository-wide rule is [Document naming standard](../DOCUMENT_NAMING.md).
The following report-specific shorthand is kept here for convenience.

Use `oneflow-<topic>-<document-type>[-YYYYMMDD].<ext>`:

- use a stable name such as `current` for actively maintained reports;
- keep the date for historical snapshots and point-in-time evidence;
- keep Markdown as the source of truth when a matching HTML report exists.

The residual baseline definitions remain under `test/baselines/`; they are test inputs,
not project reports. Raw CI/Slurm logs, environment dumps, temporary run directories,
credentials, and unredacted private cluster metadata remain outside the public repository.
