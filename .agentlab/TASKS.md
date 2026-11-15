# OneFLOW Task Board (human-maintained)

Current mode: **audit (read-only)**. Workers must not edit solver code.

## runtime-io

- A1 [todo] Read-only audit of the runtime/case-I/O boundary.
  Read the repository root `AGENTS.md`, `docs/AI_CONTEXT.md`, `docs/source/runtime_architecture.rst`,
  `codes/main/`, `codes/io/`, `tests/main/`, `tests/io/`.
  Produce `.agentlab/reports/runtime-io.md` with:
  entry points and callers; state/resource ownership; current R1 gap; relevant tests;
  shared files that cannot be owned by one worker; a proposal for disjoint future scopes;
  risks; and one concrete question for the human.

## numerics

- A2 [todo] Read-only audit of the numerics boundary.
  Read the repository root `AGENTS.md`, `docs/AI_CONTEXT.md`,
  `codes/ns/`, `codes/residual/`, `codes/grad/`, `codes/solver/`, `codes/lhs/`,
  `codes/time/`, `codes/unsteady/`, `codes/update/`, `tests/solver/`, `tests/euler/`, `tests/main/`.
  Produce `.agentlab/reports/numerics.md` with:
  which responsibilities are spatial, temporal, or shared; files that mix them;
  test coverage and missing contract tests; a proposal for disjoint future scopes;
  risks; and one concrete question for the human.

## After the audits (human)

- [ ] Review both reports and choose the first two disjoint code workstreams.
- [ ] Fill `.agentlab/contracts/` with OneFLOW-specific boundaries; do not reuse the generic CFD contract unchanged.
- [ ] Update `scope/*.txt` for code work and set `AGENTLAB_MODE=worker` in `.agentlab/target.env`.
- [ ] Register at least one real numerical/regression verification case before integration.
