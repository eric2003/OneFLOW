# Refactoring Stage Record: <short title>

Status: Planned | In progress | Partial | Verified | Blocked | Deferred
Branch: `hexin`
Base commit: `<full SHA>`
Last reviewed: YYYY-MM-DD

## 1. Problem and motivation

Describe the concrete behavior, maintenance cost, or unclear contract. Explain why it matters now.

## 2. Current implementation

- Entry point:
- Direct call path:
- State/resource ownership:
- Relevant source files:
- Relevant tests:

Separate facts confirmed from source/tests from assumptions that still need checking.

## 3. Intended contract

State observable behavior that must remain true and the behavior intentionally changing. Include compatibility constraints.

## 4. Minimal change boundary

Describe the smallest useful implementation/test seam. Record alternatives considered only when they materially affect the decision. Avoid adding an abstraction without a real consumer.

## 5. Changes made

- Source/test files changed:
- Commit(s):
- Behavior changed:
- Behavior deliberately left unchanged:

## 6. Verification evidence

| Check | Exact command or workflow | Revision | Result | Evidence |
|---|---|---|---|---|
| Focused test | `<command>` | `<SHA>` | Not run | |
| Integration/regression | `<command>` | `<SHA>` | Not run | |
| CI / platform check | `<workflow URL>` | `<SHA>` | Unknown | |

Use precise results such as Passed, Failed, Not run, Queued, In progress, or Unavailable. Do not use a green unrelated suite as proof for this contract.

## 7. Outcome and remaining gaps

- Acceptance criteria met:
- Acceptance criteria still open:
- Known limitations or risks:
- Deferred decisions and reasons:

## 8. Next action

Name one concrete next action, its expected evidence, and what condition would allow the stage to be considered complete.

## 9. Handoff summary

In 3-6 bullets, tell the next developer or AI what is established, what is not established, and where to continue. Link the relevant architecture page and `doc/AI_CONTEXT.md`.

---

Maintenance rule: keep this record focused on durable engineering knowledge. Do not record every command, intermediate thought, or routine edit. Update `doc/AI_CONTEXT.md` only when this stage materially changes the overall project status.
