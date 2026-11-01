# Verification Contract Workbench: architecture proposal

Status: **Architecture proposal; Stage 1 compatibility adapters implemented, with declared Runtime and Case contracts added as vertical slices (see Sections 13-15).**

Target branch: **hexin**

## 1. Purpose

The existing Contract Workbench is a useful numerical-knowledge browser. It already connects turbulence and solver relations to principles, references, source files, pure APIs, documentation, theory data, and plots. The next step should not turn it into a formula-only catalog. It should become a common navigation and evidence model for different kinds of engineering contracts.

A verification contract is a reviewable statement of:

1. **Intent** — what property or behavior is required, and why it matters.
2. **Scope** — which model, subsystem, workflow, case, backend, or configuration is covered.
3. **Realization** — where the behavior is implemented.
4. **Evidence** — which tests, analytical checks, cases, measurements, or reviews support the claim.
5. **Limits** — what the evidence does not establish.

The Workbench should help developers move from a claim to its implementation and evidence, and in the opposite direction from a code/test change back to the governing intent.

## 2. Design principles

- **One shared contract envelope, multiple domain-specific evidence types.** Do not force a lifecycle state machine into a curve schema, or a numerical formula into a workflow graph.
- **Evidence before presentation.** A tab is useful only when it links to a real source, test, artifact, or explicitly marked proposal.
- **Traceability is many-to-many.** One principle can govern several implementations; one test can cover several contracts; one case can exercise multiple contracts.
- **Separate claim from evidence.** A passing test is evidence for a bounded property, not proof of the whole solver or physical model.
- **Adapters own execution.** The catalog describes a contract; a typed runner or existing test system performs validation. Never duplicate production algorithms in GUI-side Python.
- **Version persisted data.** Stable IDs and schema versions matter more than the initial file format.
- **Keep the core Qt-independent.** Catalog loading, validation, evidence resolution, and report generation should be testable without a GUI.
- **Incremental migration.** Preserve the current numerical catalog and geometry scenario workflow while extracting shared concepts only where there is a demonstrated second use case.

## 3. Shared contract model

A proposed common envelope:

| Field | Meaning |
|---|---|
| schema_version | Version of the persisted contract definition |
| id | Stable, unique contract identifier |
| title | Human-readable name |
| domain | Numerical, mesh, geometry, runtime, workflow, case, I/O, build, performance, etc. |
| kind | Formula, invariant, state-transition, scenario, case-regression, workflow, interface, or benchmark |
| status | Proposed, implemented, partially verified, verified, blocked, or deprecated |
| overview | Short statement of purpose and scope |
| principle | Governing physical, mathematical, architectural, or operational reason |
| references | Literature, standards, design decisions, issue/PR, and internal documentation |
| architecture | Owning subsystem, collaborators, lifecycle/data flow, and diagrams |
| implementation | Source paths, symbols, configuration switches, and adapters |
| documentation | Canonical Sphinx/Markdown pages |
| verification | Typed evidence records and their current outcomes |
| claims | Explicitly states what the evidence proves and does not prove |
| tags | Search/filter metadata |

The envelope should not grow a large set of optional fields for every future domain. Domain-specific details belong in a typed payload selected by kind. The UI can expose common navigation tabs and enable domain-specific views only when the payload supplies them.

## 4. Contract kinds and evidence

| Kind | Typical principle / representation | Appropriate evidence |
|---|---|---|
| Numerical formula | Equations, limiting behavior, dimensional consistency | Closed-form samples, property tests, theory curves, pure-policy unit tests |
| Numerical method | Discretization, conservation, stability, order | Manufactured solutions, refinement studies, conservation budgets, benchmark cases |
| Geometry / mesh | Topological and geometric invariants | Table scenarios, mesh checks, native C++ tests, mesh artifacts |
| Runtime lifecycle | Ownership, state machine, dependency ordering | Transition tests, injected failures, production-linked integration tests, backend matrix |
| Workflow | Preconditions, dependencies, artifact contracts, failure blocking | Valid/invalid graph fixtures, runner integration tests, persisted run records |
| Case regression | Physical setup and expected behavior | Versioned input case, run manifest, residual/field/force comparisons, tolerance policy |
| I/O / configuration | Schema, parsing, rollback, compatibility | Valid/invalid fixtures, round-trip tests, production ingestion integration |
| Performance | Workload and resource assumptions | Reproducible benchmark protocol, baseline, variance, hardware/compiler metadata |

A contract may combine several evidence types. For example, a turbulence closure can link a formula curve, pure API tests, and a full CFD regression, but these remain distinct evidence levels.

## 5. Navigation model

Proposed common navigation:

- **Overview** — intent, scope, status, owner, tags, and summary.
- **Principle** — equations or invariants, physical/architectural rationale, assumptions, and limiting cases.
- **References** — papers, standards, upstream model definitions, ADRs, issues, and internal pages.
- **Architecture** — subsystem map, state/data flow, ownership, sequence/state diagrams, and boundaries.
- **Implementation** — source files, symbols, build options, runtime configuration, and adapter details.
- **Code** — relevant production excerpts and tests. Any generated sketch must be visibly illustrative and must not masquerade as source.
- **Documentation** — canonical developer/user pages and generated excerpts.
- **Verification** — linked evidence, command/configuration, result, date/revision, artifact, and limits.

Optional views are selected by contract kind: **Formula / Curve**, **Scenario**, **Case**, **Workflow Graph**, **Mesh View**, **Benchmark**, or **Run Artifacts**. These are evidence visualizers, not mandatory tabs.

The current Overview / Formula / Principle / References / Sphinx excerpt / Code sketch / Plot layout is a valid first numerical presentation. Evolve it without removing useful current behavior.

## 6. Evidence record and confidence boundaries

An evidence record should carry at least:

- stable evidence ID and evidence kind;
- linked test target/name or case/workflow ID;
- exact command and relevant configuration;
- source revision and platform/backend where applicable;
- outcome: not-run, passed, failed, blocked, or stale;
- artifact references (logs, CSV, plots, mesh snapshots, result files);
- scope and limitations.

Do not collapse these into a single green “verified” badge. A contract can have a passing unit test and still lack integration, backend, or physical validation. Display evidence by level and freshness.

Suggested evidence ladder (not a universal strict hierarchy):

1. static/schema validation;
2. unit or analytical contract test;
3. subsystem integration test;
4. production-path integration test;
5. representative numerical/case regression;
6. cross-platform/backend or multi-rank validation;
7. independent benchmark or experimental validation, where applicable.

The ladder describes increasing scope, not automatic scientific authority. A manufactured-solution test may be stronger evidence for formal order than a large engineering case.

## 7. Architecture boundaries

Proposed dependency direction:

    Contract definitions (versioned data)
                 |
                 v
    Qt-independent catalog and schema validation
                 |
                 +---- evidence resolver ---- source / tests / docs / artifacts
                 |
                 +---- typed runner adapters ---- existing test and case systems
                 |
                 v
    Presentation: Workbench tabs and domain-specific visualizers

Responsibilities:

- **Definitions** describe claims and references; they do not execute arbitrary code.
- **Catalog** validates IDs, schema versions, paths, references, and evidence declarations.
- **Evidence resolver** checks whether linked files, symbols, tests, and artifacts exist and reports stale/missing references.
- **Adapters** invoke existing test/verification infrastructure through explicit allow-listed commands and typed arguments.
- **Runner** records immutable execution metadata and artifact locations.
- **GUI** displays and edits supported metadata; it does not own solver semantics or duplicate verification logic.

Do not add a generic plugin/runner framework until at least two real domains need the same extension point. A small registry of typed adapters is preferable to dynamic imports or executable expressions in contract data.

## 8. Relationship to existing OneFLOW components

### Numerical Contract Workbench

Current implementation: tools/verification/curve_catalog.py and tools/oneflow_gui/verification_curve_dialog.py. It is currently a Python dataclass catalog with formula text, principle/background, references, pure API names, source paths, documentation paths, code sketches, proof limits, and optional CSV/plot modules.

Keep its strong feature: every numerical entry states both proves and does_not_prove. Generalize the envelope around this, not the numerical payload itself.

### Geometry Contract Workbench

Current implementation: tools/geometry_workbench and tools/oneflow_gui/geometry_contract_dialog.py. It supports editable fixture-backed input/expected tables and topology visualization. Its README correctly states that the current stub does not execute the production geometry algorithm. Preserve that distinction. A future native adapter should call production C++ or its test target rather than reimplement the resolver in Python.

### Tests and CFD regressions

GoogleTest, Python tests, grid checks, and CFD regressions remain authoritative execution systems for their respective scopes. The Workbench should link to and orchestrate them; it should not replace them with a parallel test database.

### Runtime and workflow

Runtime lifecycle contracts can link to state-transition/integration tests and runtime architecture docs. Workflow contracts can link to versioned workflow definitions, validation results, run records, and artifacts. These are natural next domains after the numerical and geometry examples.

## 9. Recommended staged implementation

### Stage 0 — inventory and semantics (no schema migration)

- Audit current numerical entries, geometry scenarios, tests, docs, and GUI tabs.
- Identify duplicate metadata and missing/broken links.
- Agree on terminology: contract, claim, evidence, scenario, case, run, artifact.
- Add no generic execution engine.

### Stage 1 — common read-only contract envelope

- Define a small versioned data model and validator in a Qt-independent module.
- Adapt one numerical entry and one geometry scenario through read-only views.
- Keep their existing source formats as the source of truth initially; avoid duplicating editable copies.
- Test schema validation, stable IDs, path containment, missing references, and proves / does_not_prove visibility.

### Stage 2 — traceability and evidence freshness

- Link exact tests, source symbols, docs, and cases.
- Detect missing/stale references.
- Record exact revision/configuration for executed evidence.
- Add a traceability report suitable for CI.

### Stage 3 — runtime and workflow contracts

- Add a lifecycle contract for process/case/accelerator ownership.
- Add a workflow contract for dependency validation and failure blocking.
- Use existing production integration tests and workflow runner; do not create duplicate engines.

### Stage 4 — case/regression evidence

- Link case definitions, immutable input snapshots, run manifests, result artifacts, tolerance policies, and comparison summaries.
- Distinguish numerical acceptance criteria from file-existence checks.
- Make results reproducible by recording executable revision, compiler, backend, and relevant options.

### Stage 5 — interactive authoring and visualization

- Add contract editing only after schema and migration behavior are stable.
- Add state diagrams, case comparisons, mesh views, and plots as type-specific visualizers.
- Preserve a CLI/headless path so validation remains usable in CI.

## 10. Acceptance criteria

The first generalization is successful only if:

- existing numerical curve entries still render and retain their principle/reference/source/test/doc links;
- existing geometry fixture scenarios remain clearly labeled as fixture-backed;
- one contract can link multiple evidence records without duplicating source/test metadata;
- missing files, tests, and invalid schema versions produce actionable diagnostics;
- the catalog and validator tests run without Qt;
- GUI presentation does not import solver implementation or run arbitrary code from data;
- evidence reports state exact revision/configuration and distinguish missing, stale, and passing evidence;
- no claim is promoted to fully verified merely because a linked unit test passes.

## 11. Industry practices to adopt

Use practices that solve concrete problems rather than copying a heavyweight process:

- **Requirements traceability:** stable requirement/contract IDs linked to implementation and verification evidence.
- **Verification and validation separation:** numerical correctness against mathematical/manufactured references is not the same as physical validation against experiments.
- **Scientific reproducibility:** record input, executable revision, compiler/backend, configuration, and output artifacts.
- **Test pyramid with production-path checks:** fast pure tests for local properties, integration tests for boundaries, and representative end-to-end regressions for system behavior.
- **Architecture decision records:** capture context, decision, alternatives, consequences, and verification status for meaningful boundary changes.
- **Schema evolution:** explicit versions, validation, migration tests, and rejection of unsupported formats.
- **Claim scoping:** every verification result says what it establishes and what it does not establish.

Useful starting references:

- NASA CFD Verification and Validation resources: https://www.grc.nasa.gov/www/wind/valid/tutorial/tutorial.html
- NASA Turbulence Modeling Resource: https://turbmodels.larc.nasa.gov/
- C++ Core Guidelines: https://isocpp.github.io/CppCoreGuidelines/CppCoreGuidelines
- JSON Schema: https://json-schema.org/

These are starting points; each contract should cite the source relevant to its algorithm, standard, or subsystem.

## 12. Open design questions for review

1. Should common contract definitions move to versioned JSON/YAML immediately, or should the first release wrap the current Python catalog and geometry JSON without migration?
2. Should definitions be authored in-repository beside code/docs, with run evidence stored outside the source tree under the user data root?
3. Should a Case be a contract kind, or should it be an evidence fixture linked to contracts? Recommended default: the latter; cases usually exercise multiple contracts.
4. Should primary grouping be domain (Numerical / Geometry / Runtime / Workflow / Case) or lifecycle (Overview through Verification)? Recommended default: domain navigation plus shared tabs.
5. Should runners launch existing CTest/GoogleTest/CFD tools directly, or initially link to commands and display externally collected results? Start with link-only/read-only evidence; add execution adapters one at a time.

This proposal remains a staged architecture guide. Stage 1 has begun with the read-only adapters described below; it does not authorize a broad rewrite.


## 13. Stage 1 implementation status (2026-11-01)

A first compatibility slice is now present on `hexin`:

- `tools/verification/contract_catalog.py` defines immutable `ContractSummary`
  and `ContractEvidence` records.
- `adapt_theory_entry` reads the existing numerical `TheoryEntry` and
  exposes its principle, references, implementation/documentation paths,
  claims, limitations, pure API names, and theory CSV as declared evidence.
- `adapt_geometry_scenario` reads an existing geometry scenario dictionary
  without changing its JSON schema or running the geometry stub / production C++.
  It explicitly labels fixture availability and linked GoogleTest references
  as not executed.
- `validate_contract_summaries` checks stable-ID uniqueness, repository path
  containment, and existence of declared implementation, documentation,
  theory-data, and linked-test source paths.
- `tools/verification/test_contract_catalog.py` adds headless tests for both
  adapters, execution-status honesty, invalid IDs, duplicate IDs, missing
  paths, path traversal, and existing repository references.

This is deliberately only a read-only compatibility model. The existing
Numerical and Geometry GUI entry points and source formats are unchanged.
The adapter does not execute tests, record results, or claim that a linked
GoogleTest has passed.

Suggested local validation from the repository root:

```powershell
python -m unittest tools.verification.test_contract_catalog
python -m unittest tools.verification.test_curve_catalog
```

The implementation has not been locally built or tested by the assistant.
The user's previously reported 1125/1125 result predates these additions and
must not be attributed to this slice. Run the two focused commands first,
then the complete project test suite.

### Stage 1 extension: combined catalog discovery

The compatibility layer now provides `load_contract_catalog(repository_root)`.
It loads every current Numerical `TheoryEntry` and every Geometry scenario
through the existing source loaders, returning immutable summaries plus
actionable diagnostics. It does not write to either source catalog.

Canonical IDs are namespaced (`numerical/<source-id>` and
`geometry/<source-id>`); `source_id` preserves the original catalog ID for
round-tripping and navigation. Geometry fixture evidence carries the
repository-relative JSON path when loaded from the catalog.

Validation now checks:
- duplicate canonical IDs;
- repository-relative path containment and file existence;
- linked GoogleTest declaration presence in its declared C++ source file.

The new Python test module is also explicitly covered by
`tools/oneflow_gui/test_test_inventory.py`, ensuring the existing
source-based Workbench test inventory discovers it. This is discovery only;
the inventory does not claim that the test has executed.

Focused validation:

```powershell
python -m unittest tools.verification.test_contract_catalog
python -m unittest tools.oneflow_gui.test_test_inventory
python -m unittest tools.verification.test_curve_catalog
```

The assistant has not executed these tests in the user's Windows environment.

## 14. First declared Runtime contract (2026-11-01)

The first cross-domain definition is now stored in
`tools/verification/contracts/runtime/process_lifecycle.json`. It describes
the process/case lifecycle and accelerator shutdown contract using the shared
navigation envelope:

- Overview and scope;
- Principle: reverse dependency teardown and retryable resource finalization;
- References to runtime architecture and shutdown-order design records;
- Architecture: `SimuImp`, `SimuContext`, and `AccelRuntime` ownership;
- Implementation source paths;
- Documentation paths;
- Typed GoogleTest evidence links;
- Explicit `proves` and `does_not_prove` claims.

`ContractSummary` now includes common `architecture` and `status` fields.
`adapt_declared_contract` accepts schema version 1 and maps a declared
cross-domain JSON definition into the same immutable read model used by the
Numerical and Geometry adapters. `load_contract_catalog` discovers JSON
definitions under `tools/verification/contracts/` without changing the
existing Numerical Python catalog or Geometry scenario files.

This is intentionally a vertical slice, not a generic schema framework:
- the contract is data, not executable code;
- GoogleTest links are checked against source declarations but remain
  `linked_not_executed`;
- `implemented` is not synonymous with `verified`;
- no runner, GUI, or arbitrary plugin execution was introduced.

The runtime definition currently claims only process-level ordering and retry
semantics. It explicitly excludes physical CUDA/HIP driver validation,
multi-rank MPI teardown, all solver-specific release paths, and CFD numerical
correctness.

Focused validation command:

```powershell
python -m unittest tools.verification.test_contract_catalog
```

The code and test definition are committed on `hexin`; the user's local test
suite has not yet validated this slice. Next work should compare the shared
model against a second genuinely different domain (workflow or case
regression) and only then decide whether any further common fields or typed
payload interfaces are justified.

## 15. First typed Case-domain contract (2026-11-01)

The next vertical slice uses a real production regression case rather than a
synthetic example:

`tools/verification/contracts/case/laminar_flat_plate_slau2.json`

It describes `test/cfd/plateuns2dslau2`, already included by
`test/suites/cpu-serial.txt`. The payload records:

- solver: `NsSolver`;
- dimension: 2;
- physical model: laminar compressible Navier-Stokes;
- convective flux: SLAU2;
- configured run length: 50 steps;
- control scripts and grid inputs;
- five output paths from `autotest/compare.txt`;
- the canonical residual database and its (10^{-15}) absolute tolerance;
- retained Blasius skin-friction and velocity reference profiles.

### 15.1 Typed payload, not a larger universal envelope

`contract_catalog.py` now defines immutable `CaseContractData` and attaches
it optionally to `ContractSummary`. The shared fields still serve navigation
and traceability. Case-specific configuration is kept in this typed payload,
rather than adding solver, mesh, and tolerance fields to every numerical,
geometry, or runtime contract.

The adapter validates required case fields and the catalog checks that suite,
case directory, input files, analytical references, and residual database
exist. It also checks suite membership, exact equality between the declared
output set and `autotest/compare.txt`, residual database schema, case
membership, residual output channels, iteration columns, row count against
`max_steps`, and agreement of the absolute tolerance. This catches stale
metadata when either the case manifest or the accepted residual history
changes.

### 15.2 What the case evidence means

This case has three distinct scientific/engineering layers:

1. **Execution and output regression.** The established `test.py` runner
   launches the solver and compares the five declared result files against
   accepted references. The ordinary text comparison tolerance is (10^{-8}).
2. **Residual history regression.** The canonical residual database checks
   all residual channels over 50 iterations with exact iteration columns,
   zero relative tolerance, absolute tolerance (10^{-15}), and rejection of
   NaN/Inf.
3. **Analytical reference material.** Blasius skin-friction and velocity
   profiles are retained under the case's `results/` directory. The current
   runner does not automatically compare these files; the Contract records
   them as declared, not executed, references.

A passing regression therefore supports reproducibility against the accepted
baseline for this configured case. It does not by itself establish independent
physical validation, grid convergence, temporal convergence, or general
accuracy across other meshes and operating conditions.

### 15.3 Why this matters for Case versus Contract

A Case is a configured experiment and an evidence fixture. It may exercise
several independent contracts: flux consistency, wall boundary conditions,
viscous transport, residual behavior, output serialization, and runtime
lifecycle. Those contracts should link to the same case evidence rather than
copying the case definition into multiple records.

The first payload is deliberately narrow. It does not add a generic case
runner, mutate suite manifests, regenerate baselines, or execute CFD during
catalog loading. Existing `test.py`, suite files, and residual database
remain authoritative.

Industry practices applied here:
- requirements traceability through stable IDs and source/test references;
- verification/validation separation (regression reproducibility is not
  independent physical validation);
- configuration management through versioned input and baseline references;
- reproducibility metadata and explicit numerical tolerance policy;
- evidence scope and limitations rather than a single verified flag.

Focused validation:

```powershell
python -m unittest tools.verification.test_contract_catalog
```

The user's previous 1126/1126 result predates this Case-domain slice. It is
implemented on `hexin`, but must remain unverified until the focused and full
test suites pass locally.

## 16. First visible Workbench integration (2026-11-01)
The Case view now provides explicit artifact navigation:
- open the declared case directory;
- open the suite manifest;
- open the residual baseline;
- open the retained analytical reference files.

These are read-only navigation actions. They resolve declared repository paths,
check containment under the repository root, and do not launch a regression or
modify accepted reference data. UI tests mock the desktop-open boundary and
assert the exact target paths, so the test does not depend on a desktop shell.


The user requested a visible cross-domain Contract Workbench rather than
leaving the catalog as JSON and documentation only. The integration is placed
in `tools/oneflow_gui`; the repository-root `gui/` remains out of scope.

The home screen now exposes one **Contract Workbench...** entry. Its domain
hub routes to four current views:

- Formula & Knowledge reuses the Verification Curves browser.
- Geometry reuses the Geometry Contract Explorer.
- Case uses a read-only domain browser over the shared Catalog and renders
  the typed `CaseContractData` in a dedicated Case Setup tab.
- Runtime uses the same common navigation without showing Case-only fields.

The shared browser presents Overview, Principle, Architecture,
Implementation/Code, Documentation/References, and Verification. This
provides consistent navigation without forcing identical payloads across
domains. Existing specialized Formula and Geometry views are preserved.

The UI is deliberately not an execution surface: it does not run CFD or
GoogleTest and does not infer PASS from a linked declaration. Catalog
diagnostics and evidence scope remain visible. Workflow and Performance
buttons are deferred until real typed definitions and evidence exist; do not
show nonfunctional placeholder buttons.

The first UI tests exercise the domain hub, Case Setup view, and Runtime
navigation. They are registered in the Python command inventory. The GUI
integration remains a thin presentation layer over the Qt-independent
Catalog, consistent with the separation between definition, validation,
execution, and presentation.


