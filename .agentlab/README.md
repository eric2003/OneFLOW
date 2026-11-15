# OneFLOW target

This target installs the agent-lab collaboration pack into a OneFLOW checkout.
It does **not** modify solver code.

## Why audit first

The starter's `spatial`/`temporal` split is a simplification. OneFLOW is a runtime/task
codebase: `codes/main`, `codes/io`, `codes/ns`, `codes/residual`, `codes/grad`,
`codes/solver`, `codes/time`, `codes/unsteady`, and `codes/update` share state and CMake
targets. Letting two workers edit code immediately would create overlapping scopes and
conflicts.

v0.1 therefore starts with two **read-only audits**:

| Workstream | Reads | Produces |
|---|---|---|
| `runtime-io` | `codes/main/`, `codes/io/`, `tests/main/`, `tests/io/`, `docs/AI_CONTEXT.md` | boundary map, R1 gap, shared files, proposed scopes, one question |
| `numerics` | `codes/ns/`, `codes/residual/`, `codes/grad/`, `codes/solver/`, `codes/lhs/`, `codes/time/`, `codes/unsteady/`, `codes/update/`, relevant tests | spatial/temporal/shared responsibility map, test gaps, proposed scopes, one question |

## Install

```bash
export AGENTLAB=/Users/luql/Downloads/agent/repo/open-research-agents
"$AGENTLAB/platform/scripts/install_target.sh" cfd oneflow /Users/luql/Downloads/agent/repo/OneFLOW
cd /Users/luql/Downloads/agent/repo/OneFLOW
git add .agentlab
git commit -m "agentlab: add OneFLOW target pack"
```

The install does not commit anything. `.agentlab` must be committed to `dev` before
workers start, because `scope_check` reads scopes from the base ref.

## Run the audits

```bash
cd /Users/luql/Downloads/agent/repo/OneFLOW
"$AGENTLAB/platform/scripts/preflight.sh" runtime-io
"$AGENTLAB/platform/scripts/run_worker.sh" runtime-io
"$AGENTLAB/platform/scripts/run_worker.sh" numerics
```

Audit gate:

```bash
GATE_MODE=audit "$AGENTLAB/platform/scripts/gate.sh" runtime-io
```

Audit workers commit only reports/RFC/repro. No build, no ctest, no code edits.

## Merge audit reports

```bash
cd /Users/luql/Downloads/agent/repo/OneFLOW
"$AGENTLAB/platform/scripts/integrate_audit.sh" runtime-io numerics
# review .agentlab/reports/*.md and .agentlab/RFC/*.md
git checkout dev
git merge --no-ff audit-integration
```

`integrate_audit.sh` refuses any change outside `.agentlab/reports`, `.agentlab/RFC`,
and `.agentlab/repro`. It does not build or run tests.

## First audit (2026-11-01)

Both audits have run and their reports are merged into `dev`:

- `.agentlab/reports/runtime-io.md`
- `.agentlab/RFC/runtime-io-2026-11-01-shared-case-runtime-target.md`
- `.agentlab/reports/numerics.md`
- `.agentlab/RFC/numerics-2026-11-01-boundary-contracts.md`

### Decided: no test-architecture refactor

The runtime-io RFC proposed a shared `oneflow_case_runtime` target so a test could
link `TeardownCase`. **Rejected for v0.1.** OneFLOW keeps its current test
architecture; `TeardownCase` stays a recorded gap. See
`platform/docs/decisions/0002-no-test-architecture-refactor.md`.

Workers may add tests only in the existing pattern (one target listing the
production sources it needs). They may not restructure the build or add a second
test entry point.

### Still open (human)

1. numerics: complete dual-time, or declare it unsupported and fail fast?
2. Which two disjoint code workstreams to open first, now that no build
   refactor is allowed.

## After the audits

1. Read `.agentlab/reports/runtime-io.md` and `.agentlab/reports/numerics.md`.
2. Approve disjoint code scopes and fill OneFLOW-specific contracts.
3. Set `AGENTLAB_MODE=worker` in `.agentlab/target.env`.
4. Verify `BUILD_CMD`/`UNIT_CMD` in `target.env` against your machine.
5. Register a real verification case before integration.

## Branch rules

- `BASE_REF=dev` for this local checkout.
- Never develop on or push to `master`.
- For upstream-compatible work, set `BASE_REF=upstream/hexin`.
- Workers do not push or merge; the human integrates and merges.
