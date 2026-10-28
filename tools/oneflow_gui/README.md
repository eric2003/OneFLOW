# OneFLOW Workbench

For the user-manual entry in the main OneFLOW documentation, see [OneFLOW Workbench User Guide](../../doc/source/workbench_user_guide_zh.rst). This README remains the detailed module-level reference for current UI behavior and service tests. The architecture overview is [Workbench Architecture](../../doc/source/workbench_architecture.rst). Case parameter file ownership, runtime database limitations, and the proposed schema-backed editor are documented in [Case Parameter Design](../../doc/source/case_parameter_design.rst).

The desktop application starts with a minimal Workbench home containing a single **Project** button. Open it to access **New Project...**, **Open Project...**, **Open Multi-Case Workspace...**, **Create Multi-Case Workspace...**, **Run History...**, and **Saved Workflows...** in a vertical list. The home page places **Tests...**, **Mesh Layout Editor...**, **Background Tasks...**, and **Configuration...** as separate vertical buttons directly below **Project**. **New Project...** creates a project from a template in a user-selected workspace directory; the initial templates are Blank project, Mesh conversion (`grid/mesh_workflow`), Mesh generation (`grid/minimal_rectangle`), and RAE2822 grid and CFD workflow. The RAE2822 template copies the mesh-conversion case and Roe-SA solver case into one workspace and includes `workflow.json`, which declares the generated-grid artifact and the solver dependency. Existing templates are copied and never edited in the source tree. The Blank project includes `.oneflow-project.json` metadata and a README. **Open Project** accepts any existing directory, including projects outside the repository and projects not registered in the built-in case catalog. The project window identifies Grid or CFD/flow projects where possible; unknown projects can still be opened. Projects inside the OneFLOW source repository must be copied to a workspace before running, and copy destinations must remain outside the source repository. User projects outside the repository can still be run in place or copied to a workspace. **Development Tests**, **Unit Tests**, and **CFD/Grid Regression** are task types with different execution strategies. The **Tests** entry groups unit tests, CFD/grid regression, and Development Tests. Development Tests is opt-in and unchecked by default. Its **Select cases...** dialog starts with all catalog cases selected, so checking Development Tests and pressing **Run** is enough; open the selector only to change the selection. Cases are copied to the fixed protected test workspace before each batch, and per-case results remain visible in the Tests dialog. The **Configuration...** dialog stores a **Run mode** preference with **Run single case**, **Multi-case run**, and **Dependency workflow** values. When **Saved Workflows...** opens, **Multi-case run** and **Dependency workflow** create a matching new draft; **Run single case** leaves the saved-workflow list unchanged. Use **Open Project...** for a single project. A versioned `.ofw` workspace file records case directories, execution mode, and dependency edges; use **Project > Open Multi-Case Workspace...** or **Project > Create Multi-Case Workspace...** to edit it. The `.ofw` editor can open the existing Saved Workflows runner through **Save and Run...**; artifact handoff still requires a full graph definition in Saved Workflows. Project discovery, catalog definitions, workflow definitions, run workspaces, and run history have distinct responsibilities.

## Workspace files

A workspace is stored as a readable JSON document with the `.ofw` extension (schema version 2). Version-1 files without artifact fields load and migrate on save. Its schema records a workspace ID and name, an execution mode (`single_case`, `sequential_cases`, `independent_cases`, or `dependency_workflow`), and an ordered list of case IDs and paths. Paths are stored relative to the workspace file when possible. Dependency workflows record `depends_on` case IDs; validation rejects unknown dependencies and cycles. The independent mode rejects declared dependencies, while sequential mode uses the case list order.

The workspace editor creates, loads, validates, edits, and runs these definitions. **Save and Run...** saves the `.ofw` file and opens an unsaved draft in Saved Workflows; execution uses the existing runner and run-history repository. Single-case, sequential, and independent modes map to one OneFLOW multi-directory process with the declared case order. Dependency mode maps IDs, directories, and edges to the existing graph runner, which stages cases and blocks dependents after failures. Dependency case rows accept `outputs` and `inputs` JSON bindings using the same fields as Saved Workflows graph definitions. Output IDs must be unique; an input must reference an output from a direct dependency. Artifact and parameter-file paths are relative to their case directory. The adapter preserves these bindings so the graph planner configures producer output settings, stages generated files, and overrides consumer settings. Version-1 definitions without `inputs` or `outputs` migrate on save; artifact fields require schema version 2.

## Launch

From the repository root, using Python 3.10 or newer:

```powershell
python -m pip install -r tools/requirements-gui.txt
python -m tools.oneflow_gui
```

Choose **Debug** or **Release** in the **Configuration...** dialog, opened from the button below **Project**; build and unit-test tasks, including Development Tests build actions, use that configuration. The executable field stays as a single compact line. Existing executables are
detected in this order: repository Debug build, repository Release build,
`C:/dev/OneFLOW/bin/OneFLOW.exe`, then conventional Windows installation
folders under `ProgramFiles`, `ProgramFiles(x86)`, and `LOCALAPPDATA`.
The first detected path is selected automatically; use the dropdown to switch
between detected copies, type a path directly, or use **Browse...** for an
unlisted location. No full-disk search is performed.

From **Open Project**, choose any existing project directory and then choose the executable and execution location. This path does not require catalog registration. To run built-in development cases, open **Tests**; Development Tests is unchecked by default. When enabled, all catalog cases are already selected; use **Select cases...** only to run a subset. The **Run mode** preference shown in **Configuration...** selects the initial draft when **Saved Workflows...** opens: multi-case creates a multi-directory draft, and dependency workflow creates a case-graph draft. It does not run the workflow automatically. Catalog cases are always copied to a run workspace before execution. The queue displays process output and per-case status. **Stop queue** terminates the active process and prevents pending cases from launching.

Catalog case execution always uses a workspace copy:
- **Copy to OneFLOW default run workspace** stages each case under
  `~/.oneflow/executions/<run-id>/<case-id>`.
- **Copy to specified run directory** stages each case under the chosen directory.
Both the source repository and its built-in case data are protected from in-place execution; workspace destinations inside the source repository are rejected.

The application invokes the standard command shape:
`OneFLOW.exe r <case-directory>`.

Use the **Tests...** button below **Project** to run the existing CTest and Python unit-test commands by default, enable the optional CFD/grid regression suite, or run the selected Development Tests cases. The shared action bar uses **Build**, **Rebuild**, **Run**, **Stop**, and **Close**. The Debug/Release selection is made in Configuration and passed into the Tests dialog. Regression runs invoke `test/test.py` and an existing manifest from `test/suites/`; the GUI does not implement a second test engine. Use **Project > Run History...** to inspect Case Workbench and workflow run records and logs, edit notes, delete individual records, or clean matching history by workflow, status, and timestamp. These history operations never delete run artifacts.

Use the **Mesh Layout Editor...** button below **Project** to open the mesh layout project editor. Choosing
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

The home screen also offers **Geometry Contracts...** for table-level interface dual-map / split scenarios (`tools/geometry_workbench`, stub backend today). This is separate from running GoogleTest.

The Workbench presents **Tests** as the entry point for unit tests, CFD/grid regression, and Development Tests. Development Tests are opt-in, and the case selector defaults to all catalog cases. `task_types.py` identifies which tasks use a user workspace; Unit Tests run from the source build tree. Test commands carry a task type, and workers publish Qt-independent `TaskRunResult` values with status, exit code when available, execution directory, and a concise message. Results are retained in run history alongside detailed case outcomes. The **Run mode** preference selects the initial draft when **Saved Workflows...** opens: multi-case creates a multi-directory draft, and dependency workflow creates a case-graph draft. It does not run the workflow automatically. Launch single projects through **Open Project...**. The multi-case workspace editor loads, edits, validates, and saves `.ofw` definitions; it does not launch runs. Use the workflow editor to create, edit, save, delete, and launch named workflow definitions.

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
- `workflow_graph_runner.py`: stages every case under a run-specific workspace, executes cases in dependency order through an injected callback, resolves declared artifact inputs against staged producer outputs, blocks consumers with missing inputs or unsuccessful dependencies, and lets independent cases continue.
- `workspace_workflow_adapter.py`: converts `.ofw` definitions into an unsaved Saved Workflows draft while keeping path resolution and execution-mode mapping independent of Qt.
- `workflow_runtime_overrides.py`: applies string settings to a parameter file in a run-specific case copy, preserving existing declaration comments.
- `workflow_case_runtime.py`: copies one case into the workflow run directory, creates configured output directories, applies producer output settings, and applies declared artifact bindings to staged consumer parameter files. The graph runner uses this service and is called by the project runner, background project runner, and Saved Workflows dialog. The dialog shares graph execution but owns its own worker and run-history updates.
- `examples/rae2822_multicase.json`: example dependency contract for a workspace whose base directory contains `rae2822grid` and `rae2822_roe_sa`. It configures the producer to write `results/rae2822_2d.ofl` and binds the solver to that exact artifact. Select `D:\\OneFLOWWorkspace\\workdir2026\\multicases` as the base directory for the user's local workspace.
- `test_workflow_case_runtime.py`: tests staged grid binding and preservation of source case files.
- `test_workflow_runtime_overrides.py`: tests existing setting replacement, missing `sourceGridCaseDir` injection, quoting, and invalid setting names.
- `test_workflow_graph_runner.py`: tests isolated staging, staged producer-consumer artifact handoff, execution order, independent branches, failure propagation, cancellation, and missing output files.
- `test_workspace_workflow_adapter.py`: tests conversion of all four `.ofw` modes, path resolution, order preservation, and dependency-edge preservation.
- `test_workflow_dependencies.py`: tests for dependency roles, ordering, missing references, and cycles.
- `test_workflow_plan.py`: tests for producer-consumer artifact binding, both `grid` and `results` paths, and invalid dependency/path declarations.
- `workflow_execution.py`: shell-free construction of one-process,
  multi-directory command arguments, preserving the user-specified order.
- `home_window.py`: Workbench activity selector with a Tests entry, configuration, project actions, and background-task management. The Tests dialog provides optional case selection for Development Tests.
- `task_types.py`: shared task categories, explicit execution-location policy, canonical run statuses, and the immutable Qt-independent `TaskRunResult` contract (task type, final status, return code, execution directory, and message).
- `test_task_types.py`: tests task classification, execution-location policy, and normalization of GUI-facing status labels to persisted states.
- `workflow_repository.py`: normalizes legacy run-status labels when loading history while continuing to persist canonical lowercase states.
- `main_window.py`: Case Workbench widgets, execution-location selection, and the sequential case queue worker.
- `mesh_layout_dialog.py`: template-based geometry, mesh-count, spacing, and boundary editor with SVG preview.
- `tools/mesh_layout/editor_document.py`: Qt-independent JSON loading, validation, compilation, and saving shared by the editor.
- `workflow_dialog.py`: workflow editor, single-process multi-directory and dependency-graph execution, and persisted run-history presentation.
- `test_services.py`: standard-library tests for the non-UI services.
- `test_workflow_execution.py`: tests for argument order, spaces, and invalid input.
- `workflow_repository.py`: versioned workflow-definition and run-history persistence, independent of Qt.
- `user_data_paths.py`: centralized cross-platform user data paths with the optional `ONEFLOW_HOME` override; the default is derived from the current user home directory.
- `application_settings.py` and `configuration_dialog.py`: persist the source root, executable, run mode, and build configuration in the per-user `settings.json`. These settings are grouped in the Workbench Configuration dialog and include reset-to-defaults support. `settings_dialog.py` is a legacy source-root-only dialog and is no longer used by the Workbench home screen.
- `test_user_data_paths.py`: tests default paths, overrides, and side-effect-free path resolution.
- `test_workflow_repository.py`: persistence tests for version checks, CRUD behavior, and history retention.

The GUI is an optional tool and is not a runtime dependency of OneFLOW. Core
case discovery and process execution remain separable from the GUI framework.
The shared `TaskRunResult` model establishes a common outcome contract without
merging the case queue and test-command workers; each runner retains its own
execution and cancellation semantics.

## Service tests

From the repository root:

```powershell
python -m unittest tools.oneflow_gui.test_task_types
python -m unittest tools.oneflow_gui.test_services
python -m unittest tools.oneflow_gui.test_user_data_paths
python -m unittest tools.oneflow_gui.test_application_settings
python -m unittest tools.oneflow_gui.test_workspace_workflow_adapter
python -m unittest tools.oneflow_gui.test_workspace_repository
python -m unittest tools.oneflow_gui.test_workspace_dialog
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

## User data directory

The Workbench keeps user-managed data under `%USERPROFILE%\.oneflow` by default
(for example, `C:\Users\eric\.oneflow`). Set `ONEFLOW_HOME` to use another
root. The source checkout remains separate from this directory.

- `settings.json`: Workbench preferences.
- `projects/`: default destination for newly created projects.
- `workspaces/`: saved multi-case `.ofw` workspace definitions.
- `workflow-data/definitions/`: saved workflow definitions.
- `workflow-data/runs/`: workflow run history.
- `workflow-data/executions/`: isolated project/workflow run copies and generated results.
- `mesh-layout-test/grid/`: fixed Development Tests execution directory, cleared as a whole and restaged with selected source cases at the start of each run batch; no per-run directory is created and source cases remain untouched.
- `mesh-layout-workspace/`: mesh layout projects, specifications, and outputs.
- `cfd-regression/`: staged CFD regression inputs.
- `background-tasks/`: persistent background task records.
- `logs/`: application logs.

The Workbench does not move or overwrite existing source projects automatically;
choose an existing project directory to open it, or create a new project under
`.oneflow/projects`.
