# OneFLOW Workbench

The desktop application starts with a Workbench activity selector. **Open Project** accepts any existing directory, including projects outside the repository and projects not registered in the built-in case catalog. The project window identifies Grid or CFD/flow projects where possible; unknown projects can still be opened. Projects inside the OneFLOW source repository must be copied to a workspace before running, and copy destinations must remain outside the source repository. User projects outside the repository can still be run in place or copied to a workspace. The **Case Database** is a separate Workbench activity, like **Run tests**. The **Run mode** control groups **Run single case**, **Multi-case run**, and **Dependency workflow**. A versioned `.ofw` workspace file records case directories, execution mode, and dependency edges; use **Open Workspace** or **Create Workspace** to edit it. Project discovery, catalog definitions, workflow definitions, run workspaces, and run history have distinct responsibilities.

## Workspace files

A workspace is stored as a readable JSON document with the `.ofw` extension. Its schema records a workspace ID and name, an execution mode (`single_case`, `sequential_cases`, `independent_cases`, or `dependency_workflow`), and an ordered list of case IDs and paths. Paths are stored relative to the workspace file when possible. Dependency workflows record `depends_on` case IDs; validation rejects unknown dependencies and cycles. The independent mode rejects declared dependencies, while sequential mode uses the case list order.

The current workspace editor creates and edits these definitions. Execution is still routed through the existing project and workflow runners; the next integration step is to make a saved workspace the single entry point for launching all four execution modes.

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

From **Open Project**, choose any existing project directory and then choose the executable and execution location. This path does not require catalog registration. To run an example from the built-in case database, open the separate **Case Database** activity; the case selector is restricted to one catalog entry. In **Run mode**, choose **Run single case** to select a user project directory, **Multi-case run** for an ordered run of independent directories, or **Dependency workflow** for cases with declared producer/consumer relationships. Catalog cases are always copied to a run workspace before execution. The queue displays process output and per-case status. **Stop queue** terminates the active process and prevents pending cases from launching.

Catalog case execution always uses a workspace copy:
- **Copy to OneFLOW default run workspace** stages each case under
  `~/.oneflow/executions/<run-id>/<case-id>`.
- **Copy to specified run directory** stages each case under the chosen directory.
Both the source repository and its built-in case data are protected from in-place execution; workspace destinations inside the source repository are rejected.

The application invokes the standard command shape:
`OneFLOW.exe r <case-directory>`.

Use **Tests...** to run the existing CTest and Python unit-test commands by default, or enable the optional CFD/grid regression suite. Regression runs invoke `test/test.py` and an existing manifest from `test/suites/`; the GUI does not implement a second test engine. Use **Run History...** to inspect Case Workbench and workflow run records and logs, edit notes, delete individual records, or clean matching history by workflow, status, and timestamp. These history operations never delete run artifacts.

Use **Mesh Layout...** to open the mesh layout project editor. Choosing
**Single rectangle**, **Single rectangle with vertical midpoints**, **Two-block turbulent flat plate**, **Two-block laminar flat plate**, or **Two-block laminar flat plate (10 edges)** copies the complete template
project once into `~/.oneflow/mesh-layout-workspace/projects`; later openings
reuse that project without overwriting edits. Each project keeps its JSON,
scripts, and generated layout together, so **Save specification** updates the
project JSON and **Compile layout** writes the layout under that project's
`grid` directory. Set `ONEFLOW_HOME` to move the workspace. Use **Load JSON...**
to open a user project at another location; that project remains in place and
is writable as long as it is outside the OneFLOW source tree. Checked-in source
files are never writable through the editor. The editor supports the single-rectangle, single-rectangle-with-vertical-midpoints, two-block turbulent flat-plate, standard two-block laminar flat-plate, and ten-edge two-block laminar flat-plate templates. The vertical-midpoint rectangle preserves four topological side positions while compiling six atomic line segments; the midpoint nodes are part of the template topology. The ten-edge laminar template splits the three vertical topological sides at y=0.02, with 41 points (40 uniform intervals) in each lower segment. The editor exposes geometry, segment point counts, spacing, and boundary types, and previews the layout as SVG. Boundary
selectors show concise names by default; hover for details or enable
**Show boundary codes** to show legacy numeric codes. The flat-plate preview
labels every edge with its edge name, boundary condition, and point count.
**Reset flat-plate boundaries** restores the standard two-block defaults.

The Workbench presents **Case Database** as a standalone activity alongside **Run tests**. The **Run mode** selector groups three execution choices: **Run single case** opens a user-selected project, **Multi-case run** opens a multi-directory workflow draft, and **Dependency workflow** opens a dependency-graph draft. Existing saved workflows remain available in the workflow editor. These are intentionally different execution semantics. Use the workflow editor to create, edit, save, and delete named workflow definitions. The
**Multi-directory** mode preserves the existing behavior: each directory is
entered on a separate line, and one OneFLOW process receives
`OneFLOW.exe r <directory-1> <directory-2> ...`. The **Dependency graph** mode
accepts a JSON array of case objects or a complete saved workflow definition.
It stages each case in an isolated run workspace, runs dependencies in order,
checks declared output files, and passes artifact paths to consumer settings.
The graph base directory contains the configured relative case directories.
Workflow definitions and run history are stored under
`~/.oneflow/workflow-data` (the user's home directory); on Windows this resolves
under the user's profile directory. Set `ONEFLOW_HOME` to use a different data
root. The existing `workflow-data` location is retained, so current definitions
and history remain available. Settings, logs, and execution workspaces have
separate paths under this root. Deleting a definition does not delete its
history. The GUI displays recent run status and messages.

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
- `home_window.py`: Workbench activity selector with standalone Case Database and Run tests entries, plus a Run mode selector for single-case, multi-case, and dependency-workflow execution.
- `main_window.py`: Case Workbench widgets, execution-location selection, and the sequential case queue worker.
- `mesh_layout_dialog.py`: template-based geometry, mesh-count, spacing, and boundary editor with SVG preview.
- `tools/mesh_layout/editor_document.py`: Qt-independent JSON loading, validation, compilation, and saving shared by the editor.
- `workflow_dialog.py`: workflow editor, single-process multi-directory execution, and persisted run-history presentation.
- `test_services.py`: standard-library tests for the non-UI services.
- `test_workflow_execution.py`: tests for argument order, spaces, and invalid input.
- `workflow_repository.py`: versioned workflow-definition and run-history persistence, independent of Qt.
- `user_data_paths.py`: centralized cross-platform user data paths with the optional `ONEFLOW_HOME` override.
- `test_user_data_paths.py`: tests default paths, overrides, and side-effect-free path resolution.
- `test_workflow_repository.py`: persistence tests for version checks, CRUD behavior, and history retention.

The GUI is an optional tool and is not a runtime dependency of OneFLOW. Core
case discovery and process execution remain separable from the GUI framework.

## Service tests

From the repository root:

```powershell
python -m unittest tools.oneflow_gui.test_services
python -m unittest tools.oneflow_gui.test_user_data_paths
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
