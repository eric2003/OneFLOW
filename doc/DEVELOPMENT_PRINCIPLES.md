# OneFLOW Development Principles

These principles are the repository's durable engineering defaults. They are intended
to prevent previously solved defects from being reintroduced during new features,
refactoring, test maintenance, or AI-assisted development.

They are defaults, not excuses to ignore a concrete requirement. If a change must
deviate, explain why, keep the exception narrow, and add evidence for the behavior.

## 1. Preserve solved knowledge

- Before changing a subsystem, inspect its existing helpers, tests, architecture notes,
  and known regression history. Do not implement a second local solution before checking
  whether the repository already has a shared one.
- A bug fix is not complete when the immediate symptom disappears. Identify the underlying
  contract and add a regression test at the narrowest useful boundary.
- When a fix reveals a reusable lesson, update this document or the owning architecture
  page in the same logical change. Record the symptom, root cause, durable rule, and
  regression test. Do not rely on chat history or a commit message as the only record.
- Prefer updating an existing authoritative document over creating overlapping notes.
  Keep chronological handoff notes concise; do not turn them into the only place where
  permanent rules live.

## 2. Keep the source tree clean and read-only during execution

- Treat checked-in source and case-definition directories as immutable inputs during
  builds, tests, GUI operations, and runtime execution.
- Put generated meshes, logs, reports, temporary files, caches, and run results in an
  explicit build or user workspace. On Windows, the normal user workspace is under
  `%USERPROFILE%\\.oneflow`; on other platforms, use an appropriate user data/cache
  location or a caller-supplied workspace.
- Use the platform's normal temporary-directory facility (for example Python
  `tempfile.TemporaryDirectory()`) for disposable data. Never move temporary data into
  the repository merely to keep it on the same drive as a source path.
- Do not assume a temporary directory and the repository are on the same filesystem or
  drive. If a tool requires same-filesystem staging, create the staging directory under
  the designated external workspace or destination filesystem, and document that
  requirement. Do not silently fall back to writing in the source tree.
- Tests must create and clean their own temporary workspaces. A test must not depend on
  stale output from a previous run.

## 3. Use portable path semantics

- Prefer cross-platform standard-library path APIs (`pathlib`, `os.path`, and
  `tempfile`) over platform-specific path manipulation.
- Do not compare paths as raw strings, use string-prefix checks for containment, or
  assume that `.resolve()` makes two aliases textually identical.
- When both paths are expected to exist and the question is whether they identify the
  same filesystem object, use `os.path.samefile(a, b)`. This handles filesystem aliases
  such as Windows 8.3 short names and symbolic links where supported.
- For normalized path policy, containment, and paths that may not yet exist, use the
  established path utility in the owning layer (currently
  `tools/oneflow_gui/path_utils.py` for Workbench/project operations). Reuse its
  `paths_equal`, `is_path_within`, and `relative_path_within` semantics instead of
  reimplementing them.
- Do not make core grid/runtime code depend on a GUI package just to reuse a helper.
  Keep utilities in the layer that owns the behavior; move one to a neutral shared
  package only when multiple independent layers genuinely need the same contract and
  the dependency direction remains valid.
- Isolate unavoidable platform-specific code behind a small, tested compatibility
  function. Do not scatter Windows API calls or OS-specific branches through business
  logic.

### Path regression to remember

A Windows CI test once compared a generated path using the normal profile path with an
expected path expressed through an 8.3 short-name alias. Both paths referred to the same
file, but direct `Path` equality (and normalization alone) did not express filesystem
identity. The test now uses `os.path.samefile`. Preserve that distinction: path spelling,
normalized path policy, containment, and filesystem identity are different questions.

## 4. Respect module ownership and dependency direction

- Trace the real entry point, callers, ownership, and cleanup path before changing an
  interface.
- Keep the C++ runtime independent of Python and the Workbench. Python may orchestrate
  the runtime, but must not become a hidden runtime dependency or duplicate numerical
  implementation.
- Keep GUI presentation separate from domain logic. Prefer Qt-independent services and
  data models where practical.
- Reuse an existing abstraction when its contract fits. Add a new shared abstraction only
  when at least two real consumers need the same stable behavior; avoid speculative
  frameworks and duplicate helpers.
- Keep platform-specific behavior at a narrow boundary and expose portable semantics to
  callers.

## 5. Keep examples and tests honest

- Each example should demonstrate one clearly named concept at an appropriate level of
  complexity. Do not label a multi-cell or multi-feature example as a single-cell case.
- Add complexity progressively. A more complicated topology, boundary arrangement, or
  physical model should be a separate example with a name that states what it exercises.
- A passing process exit code is not automatically proof that outputs, mesh topology,
  numerical values, or physical behavior are correct. Assert the actual contract.
- Keep the repository's C++ `tests/` tree C++-only. Place Python tests in their
  designated Python test area and keep the inventories and CI registration aligned.
- Separate unit, integration, regression, numerical, and platform-specific evidence.
  Report exactly what was run; do not infer a pass from an unrelated suite or an older CI
  run.

## 6. Make each change a complete, reviewable slice

- Prefer a coherent functional slice: implementation, discriminating tests, directly
  related documentation, test inventory/CI registration, and evidence as applicable.
- Do not split one behavior change into arbitrary file-by-file commits. Do not combine
  unrelated changes merely to reduce commit count.
- Preserve existing behavior unless an intentional change is documented and tested.
- Avoid hard-coded build or test file lists when the project intentionally uses a
  discovery or glob mechanism; follow the established convention for that subsystem.
- Keep source files within the repository's encoding/style constraints, including the
  C++ ASCII-source requirement.

## 7. Required pre-change and post-change review

Before editing:

1. Confirm the target branch and current HEAD; for ongoing development, use the latest
   `hexin` unless the user explicitly specifies otherwise.
2. Read this document, `AGENTS.md`, the current `doc/AI_CONTEXT.md` handoff, and the
   owning architecture/contract page relevant to the change.
3. Search for existing helpers and regression tests before adding new path, process,
   filesystem, UI, mesh, or workflow logic.
4. State the behavior contract and the evidence that will demonstrate it.

After editing:

1. Add a regression test for the root cause, not only the visible symptom.
2. Run the narrowest relevant tests first, then broader tests appropriate to the change.
3. Check that no generated artifacts or temporary files are written into the source tree.
4. Update the owning principle/architecture/experience record when the change teaches a
   reusable lesson.
5. Report the exact commit, test commands and outcomes, CI status for that exact commit,
   and anything not verified. Never claim success without evidence.

## 8. When to record an engineering lesson

Add or update a durable lesson when any of these occurs:

- a previously fixed defect is reintroduced;
- a platform-specific behavior causes a cross-platform regression;
- two modules independently implement the same policy;
- a test passed while failing to check the intended contract;
- a source-tree, workspace, process-lifecycle, or cleanup boundary is violated;
- an example's name or documentation misrepresents what it actually tests;
- a recurring review comment reveals a missing repository-wide convention.

Use this compact format in the owning document:

| Field | What to record |
|---|---|
| Symptom | What failed or was misleading |
| Root cause | The underlying semantic or architectural mistake |
| Rule | The reusable rule that prevents recurrence |
| Guard | The test, static check, or review step that enforces it |
| Scope | The modules or workflows to which it applies |

Do not record every transient command failure. Record the generalizable lesson, link the
fixing commit or test where useful, and update the rule if later evidence changes it.
