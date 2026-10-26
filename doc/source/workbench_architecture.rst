OneFLOW Workbench Architecture
==============================

Purpose and status
------------------

The Workbench is the engineering environment around OneFLOW case preparation and execution. Its purpose is broader than launching a solver process: it should make projects, run locations, workflows, execution history, mesh checks, and generated artifacts easier to inspect and reproduce.

This page is the curated architecture overview. It consolidates the bilingual engineering memos without treating every proposed capability as implemented. The code and tests for the selected `hexin` revision are the authority for current behavior.

Core domain concepts
--------------------

Project
   A named engineering context that organizes case inputs, configuration, mesh assets, outputs, and project metadata. A project may be opened from an existing directory or created from a template; templates should not be mandatory.

Case
   An existing case directory and its configuration that can be selected for execution. A case catalog helps discover and describe cases; it does not determine where a run executes.

Workspace and run directory
   A writable location for a particular execution, including its input copy, logs, and generated artifacts. The UI must state whether a run will use the source case directory, the default per-user runtime area, or a user-selected location. Source inputs should not be overwritten by default.

Workflow definition
   A saved, named execution plan. The definition is configuration and can be edited or copied independently of any execution record.

Workflow run
   One execution instance with its own ID, timestamps, executable, arguments, status, return code, logs, step outcomes, and artifact references. Editing or deleting a workflow definition must not silently erase its run history.

Artifact
   A generated or referenced file, such as a mesh, mesh summary, flow result, log, comparison, or report. Artifacts should be traceable to the run that produced them.

Execution modes
----------------

The UI and persisted model should distinguish these modes rather than hiding them behind a single ambiguous multi-case action:

* **Single-case run:** execute one selected case at a clearly displayed run location.
* **Multi-case run:** execute a user-ordered set of cases without declared data dependencies. The UI must make clear whether the mode means separate process launches or one OneFLOW process receiving multiple case-directory arguments.
* **Dependency workflow:** execute a task graph in which inputs and outputs connect steps. A dependent step must not run when a required predecessor fails or its required artifact is missing.
* **Validation and test modes:** run project checks, unit tests, case-database checks, grid/CFD regressions, or build/test commands according to their actual contracts. These are related developer operations, but should not be conflated with ordinary case execution.

A case database is an inventory and selection mechanism, not permission to execute in the source tree. Modes that can modify files should prepare an isolated runtime copy where their contract requires it.

Architectural boundaries
------------------------

* **Domain models** represent projects, cases, workflow definitions, runs, steps, artifacts, and diagnostics without depending on Qt.
* **Repositories** load, validate, and persist versioned definitions and run records. Definitions and history are separate data.
* **Execution services** resolve paths, prepare runtime copies, construct argument lists, launch processes, collect logs, handle cancellation, and persist outcomes. Avoid shell-string interpolation for process arguments.
* **Validation services** check project configuration, input files, mesh topology/geometry, and workflow prerequisites independently from presentation.
* **Presentation** displays and edits state. Window classes should not become the sole owner of persistence, process control, or validation logic.
* **Adapters** connect solver invocation, file formats, visualization tools, and report formats to stable internal contracts.

Prefer small services that can be tested without Qt. Do not add abstractions without a concrete use case; preserve a working implementation and evolve one boundary at a time.

Persistence and reproducibility
-------------------------------

Workflow definitions should use a versioned format with validation for required fields and supported versions. Unsupported versions should fail with a useful error instead of being interpreted silently.

Each run should retain, where available:

* a unique run identifier and start/end times;
* the selected executable, arguments, configuration, and resolved paths;
* the run directory and references to the input snapshot;
* process and step status, return code, and log locations;
* outputs and artifacts produced by the run;
* the OneFLOW revision or version used, when it can be determined.

The UI should show what changed in the run directory and distinguish original files from generated or modified files when that comparison is supported. Closing a details window is a presentation action, not an implicit stop request.

Background execution and shutdown
---------------------------------

Long-running tasks need an explicit lifecycle. The product design should distinguish stopping a process, closing its details view, closing the Workbench, and allowing a task to continue in the background.

If background execution is supported, a durable task record should make the task discoverable after the main window closes and after the application restarts. On restart, the Workbench should identify registered active tasks and let the user inspect or manage them before starting a conflicting new session. This behavior must be documented as implemented only when process ownership, task registration, recovery, and tests demonstrate it.

Validation and assisted development
-----------------------------------

Mesh checks, build/test operations, and structured diagnostics should be independent services with observable results. If assisted source changes are introduced, use a reviewable sequence: collect diagnostics, explain the proposed change, show a diff, obtain approval, apply it in an isolated workspace, run relevant checks, and record the outcome. Do not assume every compiler failure is an include-path problem.

Python and C++/Qt
-----------------

Use languages according to responsibility and evidence, not a predetermined rewrite plan. Python can support workflow persistence, orchestration, tests, diagnostics, and reports. Qt provides desktop interaction; C++ may be appropriate for measured visualization or native-integration needs. Both sides should share stable, versioned data contracts rather than maintaining competing project state.

Current implementation and known gaps
--------------------------------------

The implementation under ``tools/oneflow_gui`` is an optional Python/PySide6 client. ``python -m tools.oneflow_gui`` starts the home window; the C++ solver remains a separate process. The current repository contains project opening and templates, project execution in a selected source or copied location, a test command dialog, a mesh layout editor, run history, background task records, saved workflow definitions, and a separate ``.ofw`` workspace editor.

The Saved Workflows dialog supports multi-directory execution through one OneFLOW process and dependency-graph execution through ``workflow_graph_runner.py``. Project runs with a ``workflow.json`` definition also use the graph runner. The ``.ofw`` editor now offers Save and Run: ``workspace_workflow_adapter.py`` converts single, sequential, and independent modes to the existing multi-directory workflow, and converts dependency cases to the graph runner format. The dependency conversion preserves ordering, failure blocking, and artifact input/output bindings. The version-2 ``.ofw`` schema migrates version-1 definitions that do not contain artifact fields. The editor exposes per-case JSON binding columns and validates artifact IDs, direct dependencies, case-relative paths, and producer parameter settings through the canonical workflow-plan builder. It opens an unsaved Saved Workflows draft; the existing workflow dialog persists its definition and run history when the user starts execution.

``ApplicationSettings`` persists a ``run_mode`` value and ``ConfigurationDialog`` lets the user select it. When the user opens Saved Workflows, the home window maps ``multi_case`` to a new multi-directory draft and ``dependency_workflow`` to a new case-graph draft; ``single_case`` leaves the saved-workflow list unchanged. The mapping is kept in a Qt-independent helper and covered by unit tests. The preference does not yet execute a run automatically.

The Chinese user entry point is ``workbench_user_guide_zh.rst``. ``tools/oneflow_gui/README.md`` remains the detailed module reference for button behavior, data paths, service modules, and test commands. Update the user guide when visible behavior changes, the README when detailed operations or module ownership change, and this architecture page when a service boundary, persistence contract, or lifecycle changes.
Incremental implementation strategy
-----------------------------------

1. Keep domain data, workflow definitions, run history, and execution services independently testable.
2. Make run location and source-versus-copy behavior explicit in every relevant entry point.
3. Distinguish separate case launches, single-process multi-directory runs, and dependency-graph execution.
4. Add focused tests for persistence, path resolution, process lifecycle, failure handling, and artifact tracking.
5. Integrate mesh checks, case-database operations, build/test, and regression modes using their actual execution contracts.
6. Add visualization, comparison, reports, background-task recovery, and assisted diagnostics only through staged implementations with explicit tests.

This roadmap describes intended direction, not a blanket claim that all listed capabilities exist today.

Parameter editing
-----------------

See :doc:`case_parameter_design` for the proposed schema-backed GUI direction and
its migration stages. Current workflow artifact overrides remain a narrower adapter.
