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
