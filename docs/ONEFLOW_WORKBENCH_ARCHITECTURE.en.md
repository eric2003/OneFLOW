# OneFLOW Workbench Architecture Memo

## 1. Product intent

OneFLOW Workbench is an engineering workbench built for OneFLOW, not merely a GUI for launching a few cases. It should evolve alongside the solver and mesh infrastructure, gradually supporting project management, workflow orchestration, execution tracking, mesh validation, result analysis, flow visualization, data comparison, report generation, and controlled development diagnostics.

This memo records long-term intent, not a claim that every capability is implemented. Architecture should be revised as real use cases, tests, and maintenance costs provide evidence.

## 2. Core concepts

- **Project**: the boundary that organizes inputs, configuration, meshes, solver results, and project metadata. It can be opened from an existing directory or created from a template; templates are optional.
- **Case**: an existing case directory and its configuration that can be selected for execution. A case catalog helps discover and describe cases; it does not determine where a run executes.
- **Workspace**: the location for writable run copies, logs, caches, and generated artifacts. Version-controlled inputs must not be overwritten by default.
- **Run location**: the directory actually used for a particular execution. The UI should state whether it is the source case directory, the default per-user runtime area, or a user-selected directory. Modes that require isolation should prepare a copy first.
- **Run mode**: single-case execution, multi-case execution without declared data dependencies, or a dependency workflow with explicit input/output relationships.
- **Workflow**: a named, editable, saveable, copyable, and deletable execution plan. Definitions are stored separately from run records; multi-case runs and dependency workflows retain distinct semantics.
- **WorkflowRun**: a distinct record of one execution, including timestamps, executable, status, log references, step outcomes, and artifacts.
- **Artifact**: a referenced output such as a mesh summary, mesh file, flow field, comparison result, log, or report.
- **Diagnostic**: an engineering or compiler issue with source, location, severity, and suggestions.
- **Comparison**: a reproducible task for comparing runs, datasets, or baseline results.

Workflow definitions and run records must be separate. Deleting a definition must not implicitly delete history, source inputs, or run artifacts. Records should preserve, where practical, the executable path, arguments, input snapshot references, and OneFLOW version used so a run can be traced.

## 3. Architectural boundaries

1. **Domain / Models** express projects, cases, workflows, runs, steps, artifacts, and diagnostics without depending on Qt.
2. **Repositories** own versioned formats, persistence, loading, and validation. Definitions and execution history are stored separately.
3. **Execution Services** build argument lists, prepare workspaces, launch processes, collect logs, handle cancellation, and record outcomes. Commands must not be assembled through shell-string interpolation.
4. **Validation and Analysis Services** independently implement project validation, mesh topology/geometry checks, solution analysis, and data comparison.
5. **Presentation** uses Qt for display, editing, and interaction. Persistence rules, process control, and scientific logic do not belong in a giant window class.
6. **Adapters** integrate file formats, visualization backends, report formats, and existing OneFLOW interfaces without binding the core model to one GUI or format.

Core services should be unit-testable without Qt wherever practical. Avoid abstraction for its own sake: start with clear data models and small services, adding extension points when real requirements justify them.

## 4. Workflow semantics

The startup experience may present an activity-selection page for single-case execution, independent multi-case execution, dependency workflows, tests, run history, and mesh-layout tools.

Execution entry points should distinguish:

- **Single Case Run**: select one case and explicitly choose the source directory, default runtime area, or a user-selected run directory.
- **Multi-Case Run**: run multiple cases without declared data dependencies in a user-defined order. The UI must state whether this means separate process launches or one OneFLOW process receiving multiple directory arguments.
- **Dependency Workflow**: schedule steps through declared artifacts and dependency edges. A dependent step must not run when a required predecessor fails or a required artifact is missing.

Workflow steps may also cover project/input validation, builds and tests, mesh summaries and topology/geometry checks, flow analysis, data comparison, and report generation. Each step defines inputs, outputs, status, and failure behavior. Stop/continue/cleanup policy after failure should be explicit and recorded.

Workflow definitions use versioned JSON. Loading validates the version and required fields; unsupported versions fail clearly rather than being silently misinterpreted. File-based configuration and GUI editing operate on the same model, avoiding two incompatible sources of state.

## 5. History, reproducibility, and deletion

- Each execution receives a unique run ID and an independent run record.
- History is browsable by time and filterable by workflow.
- Deleting a definition is independent of retaining its history.
- Run directories and input copies should be traceable and must not overwrite source cases by default.
- Status, return code, log references, step results, and artifact references are persisted rather than living only in transient UI controls.
- Removing run artifacts is an explicit, auditable operation, not a side effect of deleting a definition.

## 6. Compiler diagnostics and assisted changes

The Workbench may gradually gain compiler diagnostics and repair assistance, but it must follow a reviewable process:

1. Collect compiler output and structured diagnostics.
2. Locate relevant source files, build targets, and dependencies.
3. Explain likely causes and propose a candidate patch.
4. Show the diff and request user approval.
5. Apply the change in an isolated workspace.
6. Rebuild and run relevant tests.
7. Record the patch, verification outcome, and failure reason.

Not every compiler error is a missing include. Causes may involve target dependencies, include paths, macros, platform differences, or incorrect module boundaries. Source files must not be silently changed without a displayed diff and approval.

## 7. Python and C++/Qt co-evolution

Python and C++/Qt may evolve together, with boundaries driven by responsibilities rather than language preference:

- Python is suitable for workflow persistence, orchestration, tests, diagnostic processing, and report generation.
- Qt provides desktop interaction; C++ components are justified for high-performance visualization, native integration, or complex interaction when needed.
- Both sides share stable, versioned data contracts rather than maintaining competing project state.
- The possibility of future C++ work should not block tested Python services today, nor should all logic be permanently forced into Python.

## 8. Incremental roadmap

1. Establish versioned workflow definitions, independent run history, and Qt-independent persistence/execution services.
2. Make single-case run locations, dependency-free multi-case execution, and dependency-graph execution explicit and persist their semantics.
3. Integrate project validation, mesh summaries, topology/geometry validation, builds, and tests.
4. Add flow visualization, result comparison, baseline analysis, and report generation.
5. Introduce structured compiler diagnostics, user-approved patches, and automatic verification.
6. Evaluate C++/Qt components based on measured performance and interaction needs instead of rewriting the whole application in advance.

Each phase should have service-level tests and a clear user-visible outcome. Do not rewrite the whole UI at once or accumulate every feature in the main window.

## 9. Immediate implementation rule

The first implementation slice establishes independent persistence for workflow definitions and run history, with tests for create, load, update, delete, version validation, and history retention. The existing directory-run modes can then be represented in the workflow model. Future capabilities in this memo are design goals, not claims of already implemented functionality; implementation status must be demonstrated by code and tests.
