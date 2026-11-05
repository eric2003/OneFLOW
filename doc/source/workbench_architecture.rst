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

The Project Window now builds a read-only Shared / CFD / Grid / Other tree from
``ParameterModel`` in control-manifest order. Effective values retain source and
declaration positions; shadowed declarations remain visible. A separate Qt-independent
``ParameterConsumerCatalog`` adds a small set of source-verified C++ read sites and
labels all other keys as not yet audited rather than unused. This curated evidence is
not whole-repository static analysis and does not prove numerical applicability.
See :doc:`case_parameter_design` for the parameter model, parser boundaries, and
consumer-evidence policy. Current workflow artifact overrides remain a narrower adapter.


Contract Workbench navigation and domain payloads
-------------------------------------------------

The shared tabs (Overview, Principle, Architecture, Implementation / Code,
Documentation / References, and Verification) describe the common engineering
contract envelope. Domain-specific content remains typed: a Case exposes Case
Setup, while Runtime lifecycle contracts expose ownership and transition details
through Architecture and Verification. Do not add empty universal tabs merely for
symmetry. Add a domain-specific view only when its typed payload and a real source of
evidence exist. In particular, a Case is currently a reproducible evidence fixture
that may exercise several contracts; it is not a replacement for the contract
itself. Workflow Graph and Run Artifacts should be introduced when workflow
contracts can link actual workflow definitions and persisted run records, not as
placeholder tabs.

A declaration or source link is not an execution result. A recorded run must retain
its revision, command, configuration, platform, outcome, and artifact references;
freshness and outcome are separate properties. A stale passing run remains a
historical pass, not current verification.


Workflow contract view
----------------------

The Contract Workbench includes a Workflow domain for declared dependency
graphs. Its typed view shows case nodes, direct dependency edges, producer
artifact paths, and consumer parameter bindings. The view is read-only and
links to the checked-in workflow definition and existing plan/runner tests.

Catalog validation checks that the typed graph agrees with the source
definition and that artifact consumers name a direct producer dependency.
Topological ordering, cycle detection, runtime staging, output checks, and
failure propagation remain owned by the existing workflow planning and
execution services. A declared graph or linked test is not an execution
record; the Workbench must not infer a successful workflow run from static
catalog validation.


Workflow execution evidence
---------------------------

The Workflow contract can be opened in the existing workflow runner through
an explicit Workbench action. The action passes the contract ID into the
runner; it does not launch a process. The developer still selects the
executable and base directory and starts the run explicitly.

After real process execution, the runner appends one ExecutionRecord per
launched step. Records contain the actual command, source case, process
working directory, timestamps, exit code, repository revision, branch,
platform/runtime metadata, explicit contract attribution, and checksummed
output/log artifacts. Dependency-blocked steps and pre-launch failures remain
in workflow history but are not represented as launched processes.

The Contract Workbench Verification view displays explicitly linked runs
separately from unlinked history and labels revision freshness independently
from the recorded outcome. A successful process exit does not override a
failed declared-output check.


Each step-level ExecutionRecord also carries the parent workflow's
`invocation_id`, allowing records to be grouped back into one graph run
without changing the stable workflow-step `test_id`. This is an additive,
optional schema-v1 field so existing records remain readable.

The design borrows the SLSA provenance separation of build definition and
run details, including invocation identity, timestamps, and output subjects
(https://slsa.dev/spec/v1.2/provenance). OneFLOW's local records are not signed
attestations and do not claim trusted-builder or complete-dependency
guarantees.


Parameter contracts in the Contract Workbench
-----------------------------------------------

The Parameter domain reuses the shared contract catalog and the existing
source-backed parameter consumer index. At this revision the index contains
32 audited runtime read-site records across 23 logical parameter keys. Each
declaration links to implementation source paths and, where supported by
evidence, related Formula & Knowledge IDs.

A Parameter contract describes the meaning and implementation of a parameter;
a Project Parameter Instance describes the effective value in one concrete
case, including its source file and override order. These are separate models.
The Contract Workbench is read-only and does not replace the Project Window
parameter editor or rewrite legacy parameter files.

The Implementation / Code preview and the Open source files action use the same
deduplicated path set: contract-declared implementation files plus audited runtime
read sites and their declared downstream behavior sources. This keeps source
 navigation consistent with the Runtime Consumers evidence tab; it does not imply
 that every linked path is executed for every Case.

Consumer coverage is deliberately conservative. A listed consumer is a
source-backed read, not proof that the parameter applies to every solver task.
An unlisted parameter is not considered unused. Declared source evidence is
not an executed test result, and a formula cross-reference does not establish
numerical order, stability, or physical validity. Future work can add typed
constraints and Case-instance links without changing this ownership boundary.


Runtime parameter-consumer evidence
-----------------------------------

parameter_consumer_catalog.py is the first source-backed consumer index. It deliberately lives outside the Qt presentation layer and stores a curated set of verified C++ reads. A record includes the parameter key, source path, consumer symbol, read expression, read timing (required / optional / lazy), applicability, and behavior/evidence notes.

Representative audited examples cover maxSteps and maxTime (outer-loop termination), time_integral and rk_stage (time-integration dispatch and stage weights), mach_ref (the derived compressible-flow reference constant), vismodel (model selection), ptTol (lazy geometric-tolerance initialization), and nmg (multigrid setup). Additional groups cover dual-time/CFL indexing, local/global time-step selection, turbulence inflow and SA production limiting, and output scheduling.

This is intentionally not presented as whole-repository usage analysis. Absence from the catalog means not yet audited, not unused. A source-level read proves a consumer exists; it does not prove the parameter is active for every task, controls a particular numerical outcome in every configuration, or has a validated physical range. Those claims require applicability tracing and separate numerical verification. The catalog is tested independently of Qt so the evidence model can later feed the Project Window and Contract Workbench.


The time-integration parameter trace is also exposed in Formula & Knowledge as `time_integral_rk_stage_weight`. Its source-backed description deliberately separates the legacy stage coefficient schedule from the narrower resident Euler backend capability gate; it does not claim a complete Runge-Kutta tableau or temporal-order verification. See `doc/source/numerical/time_integral_control.rst`.

The Project Window also displays declared Formula & Knowledge IDs for parameters that have an explicit catalog cross-reference. The first mapping connects `time_integral` and `rk_stage` to `time_integral_rk_stage_weight`; this is a navigational reference, not evidence that the temporal method has a verified order or stability region.


Parameter contract grouping and discovery
------------------------------------------

Parameter contracts carry an explicit, validated group identity, display title,
and sort order in their versioned declarations. The Parameter browser groups
entries by this metadata and offers case-insensitive search across parameter
identity, title, explanatory text, references, documentation, and implementation
paths. Search filters the existing groups rather than flattening the results.

Groups describe the engineering subject, not merely the source file that
happens to read a parameter. The initial declarations distinguish Time Integration
from Compressible Flow. Laminar and Turbulence groups should be populated only
when source-backed parameter contracts for those models are added; a parameter
must not be assigned to a physical model based on its name alone. Grouping is
navigation metadata and does not imply complete parameter coverage.


### Parameter mathematical knowledge and notation mapping

A Parameter contract may declare a typed `knowledge` payload containing
multiple formula definitions, assumptions, term-level meanings, notation
aliases, OneFLOW identifier mappings, references, and explicit caveats. The
catalog validates this payload independently of execution evidence.

The Mach-number example deliberately separates four relations: local
`M = V/a`, the ideal-gas sound-speed relation, free-stream/reference
`M_infinity = V_infinity/a_infinity`, and the exact OneFLOW initialization
expression for `const_cp`. The first three describe physical concepts under
their stated assumptions; the last is a source-derived code relation. The
Workbench's Formula & Symbols tab shows the term meanings and cross-notation
mapping, exposes the LaTeX source for copying, and can render a PNG preview
using the same renderer as Formula & Knowledge.

A symbol alias is a documented common notation, not a universal equivalence.
In particular, local Mach number and reference/free-stream Mach number must
not be conflated. A code identifier is an implementation name, not a canonical
mathematical symbol. The OneFLOW source mapping is therefore explicit and
scoped. Bibliographic references and URLs are stored with the contract, while
assumptions and unresolved physical interpretations remain visible caveats.
