# Platform rules (all agents, all domains)

1. Only modify files allowed by `.agentlab/scope/<module>.txt`. Out-of-scope changes fail `scope_check`.
2. Do not modify `.agentlab/contracts/`, `.agentlab/verification/`, platform scripts, or verification baselines/tolerances.
   Need a change -> write `.agentlab/RFC/<date>-<module>-<topic>.md`, set report status `need_decision`, and stop.
3. Completion requires `gate.sh <module>` to pass as specified for the mode, and a **commit** on `agent/<module>`.
4. Do not loosen tolerances, skip tests, delete failing tests, or mark unrun tests as passed.
5. Report honestly: `blocked` when stuck, `need_decision` when a human decision is required.
6. Ask one concrete question at a time, with a recommendation and reason.
7. Do not guess content marked TODO or [CONFIRM] in contracts; write `need_decision`.
8. Every performance or accuracy conclusion must include a repro record (see report template).
9. Do not create tasks beyond the current one, or expand your own permissions or budget.
10. Do not push, merge protected branches, or submit HPC jobs (`sbatch`/`srun`).
11. Do not restructure the target's build or test architecture. Add tests only in the
    pattern the target already uses. If a gap needs a build change, write an RFC and stop.

## Workflow

read task -> read contracts -> write tests -> implement -> gate -> commit -> write report -> stop.

Do only the first `todo` task for your module in `.agentlab/TASKS.md`.

## Commit

After the gate passes:
- `git add` only files inside your scope;
- `git commit -m "<module>: <task> <short description>"`;
- **do not push**;
- record the commit SHA in your report.

## Report status

`ready` / `blocked` / `need_decision`

`ready` means "ready for integration review"; it does not mean "numerically verified".

---

# CFD domain rules

Project: a C++ CFD solver. The **target repository** defines its C++ standard, module layout,
build system, and test runner. This domain pack does not assume a particular solver structure.

Workstreams are defined by the target's `.agentlab/TASKS.md` and `.agentlab/scope/`.
`spatial` and `temporal` are common CFD workstreams, but they are examples, not a required split.

1. Do not submit cluster jobs (`sbatch`/`srun`). Local `mpirun -np 4` is allowed; write larger commands in the report for a human.
2. Changes to numerical schemes, physical models, convergence criteria, meshes, or baselines are `need_decision`.
3. Access data and parallel primitives only through the interfaces the target declares as stable.
4. Use the target's existing test runner and labels; do not invent a parallel test framework.
5. Every accuracy/convergence/performance conclusion needs a repro record; fields are in `.agentlab/repro_fields.yaml`; state the timing method (included/excluded phases, warmup, repetitions).
6. TODO (human): naming conventions, formatter, compile flags, and the real module boundaries.

## Scope reminder

- Humans own shared foundations, contracts, verification baselines, and build-system files.
- Workers own only the paths in their `.agentlab/scope/<workstream>.txt`.
- If a change crosses workstreams, write an RFC; do not widen your own scope.

---

# --- target: oneflow ---

# OneFLOW target rules

1. Read the repository root `AGENTS.md` and `docs/AI_CONTEXT.md` before anything else.
   The root file is authoritative for branch policy, source ownership, and evidence rules.
2. Base on `BASE_REF=dev` for this local checkout. Do not touch `master`. Do not push.
   For PR-compatible work, the human may switch `BASE_REF` to `upstream/hexin`.
3. v0.1 is **audit-only**: do not edit solver code, CMake, contracts, or verification.
   Only `.agentlab/reports/<workstream>.md`, `.agentlab/RFC/<workstream>-*`, and
   `.agentlab/repro/<workstream>-*` are writable.
4. OneFLOW is a runtime/task-based codebase, not the starter's clean spatial/temporal split.
   Do not assume a boundary because the task board uses those words. Propose boundaries from code.
   The generic files under `.agentlab/contracts/` are starter templates, not OneFLOW contracts;
   see `.agentlab/contracts/README.md`.
5. Keep C++ source ASCII-only. Do not submit cluster jobs (`sbatch`/`srun`).
6. Report exact paths, callers, ownership, tests, and evidence. Do not claim a test passed unless it ran.
7. Do not restructure OneFLOW's build or test architecture (`tests/`, `test/`,
   per-target source lists, CMake targets). Add tests only in the existing pattern.
   A gap that would need a build change is a report/RFC item, not a worker task.
8. Ask one concrete question at a time, with a recommendation and reason.
9. Do not mark an audit `ready` if it only restates the task. It must produce a usable boundary map.
