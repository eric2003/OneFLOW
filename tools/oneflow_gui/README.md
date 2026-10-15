# OneFLOW Desktop Case Runner

The desktop application provides a small, extensible GUI for selecting cases
from the shared catalog and running them through the normal OneFLOW executable.

## Launch

From the repository root, using Python 3.10 or newer:

```powershell
python -m pip install -r tools/requirements-gui.txt
python -m tools.oneflow_gui
```

The executable field stays as a single compact line. Existing executables are
detected in this order: repository Debug build, repository Release build,
`C:/dev/OneFLOW/bin/OneFLOW.exe`, then conventional Windows installation
folders under `ProgramFiles`, `ProgramFiles(x86)`, and `LOCALAPPDATA`.
The first detected path is selected automatically; use the dropdown to switch
between detected copies, type a path directly, or use **Browse...** for an
unlisted location. No full-disk search is performed.

Choose an external run workspace, select one or more cases, and click
**Run selected cases**. The queue runs cases sequentially
and displays process output and per-case status. **Stop queue** terminates the
active process and prevents pending cases from launching.

Each case is copied to a unique directory under the selected workspace before
execution. Generated files therefore do not overwrite the checked-in case
inputs. The application invokes the same command shape used by the existing
workflow runner: `OneFLOW.exe r <case-directory>`.

Use **Workflows...** to create, edit, save, and delete named workflows. The
**Multi-directory** mode preserves the existing behavior: each directory is
entered on a separate line, and one OneFLOW process receives
`OneFLOW.exe r <directory-1> <directory-2> ...`. The **Dependency graph** mode
accepts a JSON array of case objects or a complete saved workflow definition.
It stages each case in an isolated run workspace, runs dependencies in order,
checks declared output files, and passes artifact paths to consumer settings.
The graph base directory contains the configured relative case directories.
Workflow definitions and run history are stored under
`~/.oneflow/workflow-data` (the user's home directory); deleting a definition
does not delete its history. The GUI displays recent run status and messages.

## Architecture direction

The long-term product direction is documented in both [Chinese](../../docs/ONEFLOW_WORKBENCH_ARCHITECTURE.zh-CN.md) and [English](../../docs/ONEFLOW_WORKBENCH_ARCHITECTURE.en.md). These memos describe the intended engineering-workbench boundaries and roadmap; they do not imply that all listed capabilities are implemented yet.

## Structure

- `models.py`: immutable case and result models plus the run-state enum.
- `case_repository.py`: catalog loading and shared validation.
- `executable_discovery.py`: ordered lookup of conventional executable paths.
- `runner.py`: workspace isolation, process launch, output streaming, and
  cancellation; it has no Qt dependency.
- `workflow_dependencies.py`: Qt-independent validation and topological ordering of case dependencies; derives independent, producer, consumer, and chained roles from graph edges.
- `workflow_plan.py`: converts a case graph into a dependency-ordered plan and binds named producer artifacts to consumer settings. Each output can declare its producer `parameter_file` and `setting`, allowing the workflow to configure the actual output path as well as the consumer input. Artifact paths are relative to the producer case, so `grid/...`, `results/...`, and other nested directories use the same contract.
- `workflow_graph_runner.py`: stages every case under a run-specific workspace, executes cases in dependency order through an injected callback, resolves artifact inputs against staged producer outputs, blocks consumers with missing inputs or unsuccessful dependencies, and lets independent cases continue.
- `workflow_runtime_overrides.py`: applies string settings to a parameter file in a run-specific case copy, preserving existing declaration comments.
- `workflow_case_runtime.py`: copies one case into the workflow run directory, creates configured output directories, applies producer output settings, and applies declared artifact bindings to staged consumer parameter files. The graph runner uses this service; the legacy GUI workflow editor is not yet connected to the graph runner.
- `examples/rae2822_multicase.json`: example dependency contract for a workspace whose base directory contains `rae2822grid` and `rae2822_roe_sa`. It configures the producer to write `results/rae2822_2d.ofl` and binds the solver to that exact artifact. Select `D:\\OneFLOWWorkspace\\workdir2026\\multicases` as the base directory for the user's local workspace.
- `test_workflow_case_runtime.py`: tests staged grid binding and preservation of source case files.
- `test_workflow_runtime_overrides.py`: tests existing setting replacement, missing `sourceGridCaseDir` injection, quoting, and invalid setting names.
- `test_workflow_graph_runner.py`: tests isolated staging, staged producer-consumer artifact handoff, execution order, independent branches, failure propagation, cancellation, and missing output files.
- `test_workflow_dependencies.py`: tests for dependency roles, ordering, missing references, and cycles.
- `test_workflow_plan.py`: tests for producer-consumer artifact binding, both `grid` and `results` paths, and invalid dependency/path declarations.
- `workflow_execution.py`: shell-free construction of one-process,
  multi-directory command arguments, preserving the user-specified order.
- `main_window.py`: Qt widgets and the sequential queue worker.
- `workflow_dialog.py`: workflow editor, single-process multi-directory execution, and persisted run-history presentation.
- `test_services.py`: standard-library tests for the non-UI services.
- `test_workflow_execution.py`: tests for argument order, spaces, and invalid input.
- `workflow_repository.py`: versioned workflow-definition and run-history persistence, independent of Qt.
- `test_workflow_repository.py`: persistence tests for version checks, CRUD behavior, and history retention.

The GUI is an optional tool and is not a runtime dependency of OneFLOW. Core
case discovery and process execution remain separable from the GUI framework.

## Service tests

From the repository root:

```powershell
python -m unittest tools.oneflow_gui.test_services
python -m unittest tools.oneflow_gui.test_workflow_repository
python -m unittest tools.oneflow_gui.test_workflow_execution
python -m unittest tools.oneflow_gui.test_workflow_dependencies
python -m unittest tools.oneflow_gui.test_workflow_plan
python -m unittest tools.oneflow_gui.test_workflow_graph_runner
python -m unittest tools.oneflow_gui.test_workflow_runtime_overrides
python -m unittest tools.oneflow_gui.test_workflow_case_runtime
python tools/test_mesh_cases.py
```

The first version reports process exit status and logs. It does not yet claim
that a successful process exit proves mesh geometry is valid; result validators
will be integrated as a separate workflow stage.
