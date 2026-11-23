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
