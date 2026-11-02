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

The Tests workbench reports two separate accounting levels. **Command counts**
(Passed / Failed / Not run) describe launched command groups; **case counts**
(Passed / Failed / Skipped) describe framework-reported test cases. A command
that was never launched is not a skipped test case. Python unittest and CTest
counts remain separate on the final status line so different frameworks are
not silently conflated. CTest skipped cases are counted from its explicit
"tests did not run" section; absent that section, the skipped count is zero.


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

## 17. Evidence reference path validation (2026-11-02)

The shared catalog validator now checks repository-file evidence declared as
`output_manifest`, in addition to the existing theory-data and geometry-fixture
references. This closes a traceability gap: a case contract could previously
point its output-manifest evidence at a missing or misspelled file while the
case-specific validator independently checked the conventional
`autotest/compare.txt` path.

The contract's declared evidence reference is now validated for repository-root
containment and file existence. The existing case-specific comparison still
checks that the manifest's output set agrees with the typed case payload. These
are complementary checks:

- evidence-path validation asks whether the exact artifact named by the
  contract exists and stays inside the repository;
- case-manifest validation asks whether the artifact's declared outputs agree
  with the case contract.

A focused negative test mutates the manifest reference to a missing file and
requires an actionable catalog diagnostic. Catalog loading remains read-only;
no solver, test runner, or GUI dependency is introduced.

Focused validation:

```powershell
python -m unittest tools.verification.test_contract_catalog
```

The user reported 274/274 Python unittest cases passing before this change.
That result is the prior green baseline, not validation of this new test.

## 18. Typed evidence kind validation (2026-11-02)

Versioned contract adaptation now validates the evidence `kind` against an
explicit schema-v1 vocabulary, just as it already validates evidence status.
Previously, any non-empty kind string could enter the immutable summary. An
unknown kind would then be ignored by kind-specific validation branches, which
could make a malformed declaration appear structurally acceptable.

The current vocabulary covers numerical API/theory data, geometry
fixture/linked-test evidence, GoogleTest links, and Case suite/output-manifest/
residual-baseline/analytical-reference artifacts. New kinds must be introduced
intentionally with parser and validation behavior, documentation, and tests.

Kind and status are separate dimensions:
- kind identifies what artifact or relationship is declared;
- status describes whether it is merely declared, linked, recorded, passed,
  failed, blocked, or stale.

Recognizing a kind does not execute it or upgrade its status. This keeps the
catalog a read-only definition and validation layer.

Focused validation:

```powershell
python -m unittest tools.verification.test_contract_catalog
```

A negative test changes a valid Runtime evidence kind to an unknown value and
requires adaptation to reject it with an actionable error.


## 19. Cross-domain artifact evidence path validation (2026-11-02)

Artifact references now receive the same repository-root containment and
file-existence checks across domains. The shared path-evidence vocabulary
includes theory data, case-suite manifests, output manifests, residual
baselines, and analytical references. Geometry fixture references are checked
when they are path-shaped; an opaque scenario ID remains an identifier rather
than being mistaken for a filename.

The path list is deduplicated before validation so a Case payload and its
typed evidence can refer to the same baseline or analytical file without
duplicating diagnostics. Case-specific semantic checks remain in place:
existence is not enough to establish suite membership, output agreement,
residual schema, row count, or tolerance consistency.

Negative tests attach missing analytical-reference and residual-baseline
artifacts to a Runtime contract and require actionable diagnostics. A separate
test protects opaque geometry fixture IDs from accidental path interpretation.

Focused validation:

```powershell
python -m unittest tools.verification.test_contract_catalog
```

The validator remains read-only and Qt-independent. It reports broken
declarations; it does not execute evidence or promote a contract to verified.


## 20. Typed GoogleTest declaration matching (2026-11-02)

The linked GoogleTest source resolver now recognizes both ordinary and typed
test declaration macros: `TEST`, `TEST_F`, `TEST_P`, `TYPED_TEST`, and
`TYPED_TEST_P`. Previously the matcher accepted only the ordinary TEST
family, so a valid typed-test link could be reported as stale.

A focused test writes both typed declaration forms into a temporary C++
source file and resolves contract evidence against that file. This avoids
adding artificial production test declarations to the repository.

The boundary remains intentionally narrow: matching a source macro establishes
that the named declaration is present in the declared source file. It does not
establish test-binary registration, successful compilation, test discovery at
runtime, or execution. Those require build/test inventory and execution
records, which remain separate evidence layers.

Focused validation:

```powershell
python -m unittest tools.verification.test_contract_catalog
```


## 21. Ignore non-code text when resolving GoogleTest declarations (2026-11-02)

The linked-test source check now masks C++ line comments, block comments, ordinary string/character literals, and raw string literals before matching GoogleTest declaration macros. A declaration-shaped example in documentation text or a literal can no longer make a stale contract link appear resolved. New negative tests cover all four source forms.

This is a deliberately bounded lexical check rather than a C++ parser. It establishes only that a declaration-shaped macro invocation appears in code text in the declared source file. Macro expansion, parameterized/typed suite registration, binary discovery, compilation, and execution remain outside this resolver and must be represented by later build/execution evidence.

Focused validation:

```powershell
python -m unittest tools.verification.test_contract_catalog
```



## 22. Separate evidence reference resolution from execution status (2026-11-02)

The Qt-independent catalog now exposes `EvidenceResolution` and `resolve_contract_evidence(repository_root, summary)`. Every declared evidence item receives a static resolution state: `declared` when no supported artifact resolver applies, `resolved` when its repository file or linked GoogleTest declaration is found, and `stale` when that supported reference is invalid, missing, unreadable, or absent from the declared source.

The original `ContractEvidence.status` is retained unchanged. For example, a source declaration can be `resolved` while its declaration status remains `linked_not_executed`. This deliberately avoids manufacturing a PASS result from static source inspection. The Case-specific semantic checks remain in the catalog validator; this first resolution API does not claim to validate scientific contents of analytical data or output artifacts.

The Runtime Verification view now displays resolution status, declaration status, scope, resolved file, and any diagnostic separately. Linked-test scopes may append explanatory text after a semicolon; the resolver extracts and validates the leading C++ source path instead of treating the full annotated scope as a filename. Focused tests cover resolved, stale, and declaration-only evidence, verify declaration status is preserved, and assert the GUI presents both dimensions.

Focused validation:

```powershell
python -m unittest tools.verification.test_contract_catalog
python -m unittest tools.oneflow_gui.test_contract_workbench_dialog
```

## 23. Static CMake target linkage for GoogleTest evidence (2026-11-02)

Added a Qt-independent CMake test catalog that maps a resolved GoogleTest
source path to direct `add_executable` source declarations in
`tests/**/CMakeLists.txt`. It separately records whether the matched target
has a `gtest_discover_tests` declaration and, when present, its
`DISCOVERY_MODE`. The Runtime Verification view now displays the target,
owning CMake file, target-resolution state, and discovery declaration.

The state boundary is explicit:

- source declaration resolved: the named GoogleTest macro was found in source;
- target resolved: a direct CMake executable declaration lists that source;
- CTest discovery declared: the target has a `gtest_discover_tests` declaration;
- built, discovered at configure/test time, and executed: **not established by
  this static resolver** and must come from actual build/CTest records.

The parser deliberately does not infer target membership from variable-based
source lists, `target_sources()`, helper functions, or generated CMake. Such
cases remain unresolved until a configured build-system export or CMake File
API integration can provide authoritative target metadata. This conservative
boundary is preferable to a plausible but false green state.

Tests cover positive target/discovery mapping, a target without discovery
registration, variable-indirect source lists, unrelated files, and paths
outside the repository. The GUI integration test checks the actual Runtime
contract target and its PRE_TEST discovery declaration.

Focused validation:

```powershell
python -m unittest tools.verification.test_cmake_test_catalog
python -m unittest tools.oneflow_gui.test_contract_workbench_dialog
```


## 24. Structured execution diagnostics (2026-11-02)

The Tests command worker now turns failed Python unittest and GoogleTest
output into immutable `TestDiagnostic` records. The first diagnostic slice
captures:

- failed task label and test framework/identity;
- source file and line when emitted by the framework;
- exception type and message when present;
- selected CMake configuration, resolved working directory, UTC start and
  finish timestamps, and process exit code;
- explicit contract IDs and a conservative verification-impact field.

The Tests dialog renders the record as structured fields while retaining the
raw output and command summary. Missing information is represented as unknown.
The parser does not infer contract identity from test names or source proximity.
The current runner supplies no automatic contract mapping, so the default
impact is `unlinked`. If a future resolver supplies explicit contract IDs,
the diagnostic is marked `requires_review`; a failed test is not automatically
equivalent to disproving every claim in a contract.

This slice intentionally stops before persistent execution storage and
Verification Results aggregation. The next integration should add a versioned
execution-record repository under the user data root, preserve command and
per-test records, and resolve test links through stable evidence IDs. The
Contract Workbench should then display the latest matching execution evidence
without mutating contract declarations. Verification conclusions must be
computed from explicit contract scope, evidence kind, run status, and freshness;
a missing run remains unknown, not failed or passed.

Focused validation:

```powershell
python -m unittest tools.oneflow_gui.test_test_diagnostics
python -m unittest tools.oneflow_gui.test_test_runner
python -m unittest tools.oneflow_gui.test_test_inventory
```

The assistant has not run these commands in the user's Windows environment.


## 25. Make failed-case summaries actionable (2026-11-02)

The test worker now accumulates structured diagnostics from each failed test
command and appends a final `Failed test cases` section after the aggregate
command counts. Each listed case includes its task label, source location,
and concise failure cause. The original streamed framework output remains
above the final summary and retains the full traceback/assertion detail.
When structured records exist, the worker suppresses the older generic
unittest name list to avoid duplicate entries.

Diagnostic working-directory values are serialized with forward slashes so
the same resolved path has a stable display representation on Windows and
POSIX. The process still executes with the original native `Path`.

A regression test now verifies the final case list includes the test identity,
file/line, and assertion message. Another assertion compares normalized path
representation rather than requiring a platform-specific separator.

Focused validation:

```powershell
python -m unittest tools.oneflow_gui.test_test_diagnostics
python -m unittest tools.oneflow_gui.test_test_runner
python -m unittest tools.oneflow_gui.test_test_inventory
```


## 26. Cross-platform execution assumptions in the test workbench

The test command layer keeps argument vectors as the source of truth and uses
host-appropriate quoting only for human-readable logs. Windows-specific
executable suffixes and multi-configuration output paths are selected only on
Windows. POSIX single-configuration builds use the conventional unsuffixed
executable under `build/bin`; POSIX multi-configuration generators are also
recognized through their CMake cache metadata. The configure command no longer
injects a local Windows install prefix, leaving CMake's platform default intact.

This is a bounded compatibility improvement, not a claim of full Linux
validation. Platform branches are unit-testable through an explicit platform
parameter, while CI or a developer Linux machine must still perform a real
configure, build, CTest, and regression run to establish runtime compatibility.

## 27. Cross-platform process-tree cancellation

The GUI test worker, Case Database runner, and workflow process runner now
start each command in an isolated process group via the Qt-independent
`tools.oneflow_gui.process_control` module. Stop terminates the whole process
tree rather than only the direct child:

- POSIX: `start_new_session=True` and process-group termination with
  `os.killpg`;
- Windows: `CREATE_NEW_PROCESS_GROUP` and
  `taskkill /PID <pid> /T /F`, invoked without a shell.

This closes an execution-lifecycle gap in CFD regression. The regression
Python script may start a shell command, which may start MPI and solver
processes. Parent-only termination is not a sufficient cancellation contract.

The module has focused tests for host-specific process creation options,
Windows tree termination and fallback, and completed-process handling. Native
CI exercises a real parent/descendant process tree on both POSIX and Windows;
the Windows path uses `taskkill /T`. The new test is part of the Python test
inventory and the Ubuntu/Windows Python tooling CI matrix.

The implementation deliberately keeps platform-specific process mechanics at
the execution boundary. Solver algorithms, regression command definitions,
and task result semantics remain independent of the host OS. MPI launcher
behavior still requires an end-to-end stop test on each supported host.


## 28. Link Python process tests to the shared contract catalog

The Runtime domain now has a separate
`runtime/gui-process-execution` contract for operating-system processes
launched by the desktop workbench. It is intentionally distinct from
`runtime/process-lifecycle`, which describes accelerator backend ownership
and finalization inside the solver. Similar words such as *process* and
*lifecycle* must not collapse different resource owners into one contract.

The contract links the shared process-control shutdown policy and the three
GUI execution paths to concrete Python unittest methods. The catalog recognizes
the new `python_test` evidence kind by parsing source with Python's AST and
checking the declared test class/method. This is static evidence resolution:
it does not import or run the module, and it preserves the explicit
`linked_not_executed` status. A missing or renamed method becomes stale rather
than being silently matched by a similar test name.

The catalog tests cover both successful Python test declaration resolution and
a missing method. The contract itself is loaded and validated as part of the
combined catalog test, and that test is now included in the Ubuntu/Windows
Python tooling CI matrix. This makes the path from runtime claim to source test
navigable without conflating declaration with execution evidence.

Focused validation:

```powershell
python -m unittest tools.verification.test_contract_catalog
python -m unittest tools.oneflow_gui.test_test_inventory
```

The next stage is a versioned execution-record model. It should capture the
actual command, host/platform, revision, timestamps, exit status, diagnostics,
and artifact references, then match those records to these declared evidence
IDs. Until that exists, the catalog correctly reports the new links as
unexecuted regardless of local or CI results.

## 29. Durable execution-record foundation

The first execution-record vertical slice is implemented in
`tools/verification/execution_records.py`. It is intentionally independent of
Qt and does not replace `TaskRunResult` or `WorkflowRunRecord`:

- task results describe GUI task scheduling and completion;
- workflow records describe workflow runs and their editable history;
- verification execution records describe one test invocation and the
  provenance required to assess its evidence.

Schema version 1 captures stable test identity, framework, source path and
optional line, exact command, repository/commit/branch, configuration,
platform/runtime version, working directory, timezone-aware timestamps,
exit code, outcome, explicit contract IDs, diagnostics, and artifact paths
with optional SHA-256 checksums. Optional environment data is allowed to be
unknown; no missing field is filled by inference. A record can be current,
stale, or unknown relative to a selected commit.

The JSONL store is append-only and strict on read: malformed lines report
their line number instead of silently losing evidence. Its default path is
reserved under the user data root, separate from workflow history.

The first real execution adapter is now connected to
`TestCommandWorker` for Python `unittest` commands. It persists a record only
after a process has launched and returned, captures the exact argument vector,
timestamps, exit code, working directory, available Git/runtime metadata, and
a SHA-256-addressed stdout artifact. Pre-launch cancellation and process-start
failure do not create execution records; post-launch cancellation is recorded
as cancelled. Contract IDs remain empty until an explicit mapping is authored.

This proves the persistence boundary for one framework, but the Contract
Workbench does not yet render execution history. The next vertical slice should
display recent records and freshness in the Workbench, while keeping static
resolution separate from actual execution. Only after that path is proven
should CTest/GoogleTest and CFD regression adapters be added. A declaration,
static source resolution, build registration, execution, and passing outcome
must remain separate states.


## 16. Recorded execution history presentation

The domain Verification tab now displays immutable execution records from the
user-data JSONL store. The presentation has two deliberately separate groups:
records explicitly linked through `contract_ids`, and recent unlinked records.
Unlinked records are visible for discovery but are not attributed to the selected
contract. Each record shows outcome, commit freshness, exact command, completion
time, and artifact references.

Freshness is compared against the current checkout's Git revision. If either
revision is unavailable, the status remains unknown. This is a read-only
presentation layer: it does not infer links from test names, source paths, or
matching commands, and it does not claim that a passing run verifies the entire
contract. The integration test exercises both linked and unlinked histories and
current/stale revision states.


## 17. Explicit links from Python execution records

The Python unittest capture adapter now resolves contract IDs through the shared
catalog. It links a record only when the command resolves to one exact Python
module source file and a contract explicitly declares Python test evidence from
that same file. Directory-wide discovery, unknown files, and catalog failures
remain unlinked. The resolver is Qt-independent and has focused catalog and
capture tests. This is an evidence association, not a claim that every test in
the file was relevant to every contract or that the contract is fully verified.


## 30. Shared Workbench page composition across engineering domains

The common Workbench should standardize navigation and evidence semantics,
not force every engineering topic into a formula-shaped record or view.

### Shared page sections

The proposed shared navigation vocabulary is:

1. **Overview** - intent, scope, domain, kind, maturity, and concise summary.
2. **Principle** - governing equations, invariants, ownership rules, state
   transitions, numerical assumptions, or operational constraints.
3. **References** - papers, standards, upstream definitions, ADRs, issues, and
   internal design documents.
4. **Architecture** - module boundaries, ownership, dependency direction,
   sequence/state/dataflow diagrams, and collaboration maps.
5. **Implementation** - production entry points, configuration, adapters,
   data structures, and relevant symbols.
6. **Code** - source and test navigation with clearly separated production
   code, test code, and illustrative snippets.
7. **Documentation** - canonical developer/user pages and generated excerpts.
8. **Verification** - declared tests, actual runs, evidence artifacts, outcomes,
   freshness, and explicit limitations.

The navigation labels are shared, but a section is shown only when the contract
has meaningful content for it. Missing evidence must be shown as missing; an
inapplicable view should not become an empty placeholder.

### Domain-specific visualizers

A domain may provide one or more specialized visualizers without changing the
shared envelope:

| Domain / kind | Specialized visualizers |
|---|---|
| Numerical formula or method | Formula, derivation, theory curve, convergence plot |
| Mesh or geometry | Topology, orientation, boundary map, mesh snapshot |
| Runtime lifecycle | State machine, sequence diagram, ownership graph |
| Case / CFD regression | Case setup, input inventory, run artifacts, residual and field comparisons |
| Workflow | Dependency graph, artifact flow, failure propagation |
| Performance | Benchmark protocol, scaling plot, profile and resource summary |

A specialized visualizer is a projection of domain-owned data. It is not a new
source of truth and must not duplicate the solver algorithm in GUI code.

### Three independent status dimensions

The page must distinguish:

* **Intent status** - proposed, accepted, deprecated, or superseded design.
* **Implementation status** - absent, partial, or implemented production behavior.
* **Evidence status** - not linked, declared, resolved, executed, passed, failed,
  stale, or blocked, with the scope and revision of the evidence.

These dimensions must not collapse into one verified flag. For example, an
implemented runtime lifecycle may have passing unit tests while still lacking
multi-rank MPI or physical accelerator evidence.

### Layout and interaction proposal

The initial layout can retain the current Formula Knowledge experience as one
specialization. A shared page should place the contract identity and status
summary at the top, shared section navigation below it, and domain visualizers
inside the relevant section. Architecture diagrams should link to their source
modules and production entry points; implementation references should link to
source files and symbols; verification entries should navigate to the exact
test declaration and, when available, the immutable run record and artifacts.

For Runtime, the Architecture view may show process and case scopes as a
sequence/state diagram. For Case, the Principle view may explain governing
equations and boundary conditions while the Case view shows actual inputs and
run outputs. For Mesh, the Principle view may state orientation and incidence
invariants while the Topology view visualizes nodes, faces, cells, and boundary
ownership. This keeps a consistent mental model without pretending all domains
are numerical formulas.

### Migration rule

Do not replace the existing Numerical or Geometry dialogs in one rewrite.
Keep their current catalogs as authoritative inputs, continue using the
read-only adapters, and introduce the shared shell only after the same
navigation and evidence behavior has been exercised by at least two different
domains. Runtime and Case are the next useful examples. The common schema should
remain small; domain-specific payloads and visualizers should be typed and
validated independently.

### Acceptance checks

* Numerical Formula/Curve behavior remains unchanged.
* Runtime pages can render architecture and lifecycle diagrams without a
  formula field.
* Case pages can show physical setup, input assets, run evidence, and residual
  comparisons without requiring a formula plot.
* A section with no applicable content is omitted or clearly marked
  inapplicable, never silently treated as verified.
* Each displayed implementation or evidence link resolves to its declared
  source; execution status is never inferred from a declaration.
* The Qt-independent catalog and validation remain usable in CI.

## 31. Selective source navigation

Source-opening actions must never launch every associated file without asking.
The shared file-opening path first resolves and deduplicates repository-relative
files, rejects paths outside the repository, and excludes missing files. It then
presents a selection dialog with extended selection, Select all, and
Clear selection actions. A user can open one file, several selected files,
or all files; Cancel opens none. The same selection behavior applies to
implementation sources, documentation and references, and declared test
sources. This keeps the action predictable when a contract links multiple
production modules or tests.

The Runtime contract runtime/gui-process-execution must expose its declared
implementation paths in the Implementation / Code view and use the same
selection action as every other contract. A regression test covers visibility
of this action, the declared process-control source paths, selection cardinality,
and opening only the selected file.



## 32. Put source-opening actions in one visible action row

The Implementation / Code source action is now part of the shared action row
below the tabs, rather than being placed inside the Implementation page where
it can be missed. The row order is intentional: **Open source files** first,
then **Open documentation**, then **Open declared test sources**, followed by
Case-specific actions when the Case domain is active. The source action always
opens the implementation paths declared by the selected contract; it does not
substitute the smaller set of test-evidence paths.

The file picker still supports one, multiple, or all selected files. Source
files are launched by an explicitly selected editor process rather than by the
Windows file association for ``.py`` files; this avoids accidentally starting
the Python interpreter and closing immediately. Configuration persists a
source-editor preference: Visual Studio Code (default), Visual Studio 2026,
Visual Studio 2022, Notepad++, or the system default application. Missing
executables produce an actionable warning. Directory navigation continues to
use the operating system file manager.


## 2026-11-02 continuation: repository-case setup-only evidence

The R1 production integration target now declares a setup-only smoke test for
the checked-in test/grid/rae2822convert control scripts. The test copies the
script directory to a temporary case, uses the production Context setup and
teardown methods, asserts typed database values and generated diagnostics,
and does not resolve or execute the Grid task. This avoids source-tree output
pollution and keeps control ingestion separate from grid conversion.

The case's grid.txt and part.txt retain legacy asset names that do not all
match the case's checked-in files. Consequently, successful setup is not
evidence that the Grid task or conversion output is valid. A separate task
execution/output contract is required for that claim.

The runtime contract declaration links the exact GoogleTest name as
linked_not_executed. Static catalog resolution, CMake target declaration,
CTest discovery, actual execution, and pass outcome remain distinct evidence
states. The assistant has not run the user's Windows build or CTest.


## 25. First typed Workflow-domain contract (2026-11-02)

The first Workflow-domain slice uses the checked-in RAE2822 grid-to-solver
example:

`tools/oneflow_gui/examples/rae2822_multicase.json`

The declared contract
`tools/verification/contracts/workflow/rae2822_grid_to_solver.json`
captures its two case nodes, the direct producer-consumer dependency, the
stable `rae2822_grid` artifact ID, the producer output path, and the
consumer's `gridFileName` binding in `script/cfd.txt`.

### Typed graph payload and consistency

The shared immutable `ContractSummary` now optionally carries
`WorkflowContractData`, composed of typed step, input-binding, and
output-artifact records. This remains a domain-specific extension; no
workflow-only fields were added to numerical, geometry, case, or runtime
contracts.

The catalog checks:
- unique workflow step IDs and output artifact IDs;
- non-empty, root-relative case directories, artifact paths, and parameter
  file paths without parent traversal;
- every consumer input resolves to a declared producer artifact;
- every consumer declares a direct dependency on that producer;
- the typed graph payload matches the checked-in source workflow definition.

The catalog does not reimplement topological sorting, cycle detection, runtime
staging, process execution, or artifact existence checks. Those remain owned
by the existing workflow plan and graph runner. Their Python tests are linked
as `linked_not_executed` evidence until an actual execution record is attached.

### Workbench presentation

The domain hub now exposes Workflow. Its Workflow Graph tab presents the
producer and consumer steps, dependency edge, artifact path, and parameter
binding. A read-only action opens the exact workflow definition; the existing
Verification tab links the plan and runner tests. This is a presentation of
declared intent, not a claim that the example has been executed.

Focused validation:

```powershell
python -m unittest tools.verification.test_contract_catalog
python -m unittest tools.oneflow_gui.test_contract_workbench_dialog
python -m unittest tools.oneflow_gui.test_workflow_plan
python -m unittest tools.oneflow_gui.test_workflow_graph_runner
python -m unittest tools.oneflow_gui.test_test_inventory
```

The test inventory now asserts discovery of the new catalog and Workbench
tests. Existing CI test commands continue to execute the inventory and
registered Python test groups; no duplicate workflow-specific CI pipeline
was added. The user's 1156/1156 result is the pre-change baseline; this slice
remains unverified until the focused and full local suites pass.


## 26. Workflow execution evidence reaches the Workbench

The Workflow contract now has an explicit **Run in workflow runner...** action.
It opens the existing workflow editor/runner with the checked-in definition
and passes the contract ID as an explicit attribution. The runner still
requires the developer to choose the executable and graph base directory and
to press Run; opening a contract does not execute a case.

When a process is actually launched, the workflow completion path appends
one immutable schema-v1 ExecutionRecord per launched workflow step. Each
record retains the exact command arguments, source case directory, process
working directory, start/finish timestamps, exit code, repository revision,
branch, platform, Python runtime, recognized build configuration, and
checksummed output/log artifacts. The explicit contract ID is copied into the
record only when the run was launched from that contract's Workbench action.

Dependency-blocked steps and failures before process launch remain in the
existing WorkflowRunRecord history. They are not fabricated as process
executions. A process that exits successfully but fails the declared-output
check is recorded as failed even when its process exit code is zero. Runs
launched from the general workflow manager without an explicit contract link
remain unlinked in the Contract Workbench.

The existing Verification view displays linked and unlinked execution history
separately and compares the recorded commit with the checked-out HEAD. This
connects a real workflow run to the declared contract without introducing a
second workflow engine.

Focused validation:

```powershell
python -m unittest tools.oneflow_gui.test_workflow_execution_capture
python -m unittest tools.oneflow_gui.test_contract_workbench_dialog
python -m unittest tools.oneflow_gui.test_workflow_repository
python -m unittest tools.oneflow_gui.test_workflow_graph_runner
python -m unittest tools.oneflow_gui.test_test_inventory
```

The new test module is registered in both the explicit Python unittest registry and Python Tooling CI. The user's 1156/1156
result is the pre-change baseline; this feature slice requires its own focused
and full-suite validation.


### Execution identity and provenance boundary

Each launched step record now carries the parent workflow's stable run ID as
`invocation_id`. This makes the per-step records from one graph execution
correlatable without overloading `test_id`, which remains the stable identity
of the workflow step. Existing schema-v1 records without this optional field
remain readable.

This follows the SLSA provenance model's separation of build definition and
run details, including invocation identity, timestamps, and output subjects:
https://slsa.dev/spec/v1.2/provenance

OneFLOW's local ExecutionRecord is a developer-facing execution log, not a
signed SLSA attestation. It does not claim a trusted builder identity,
complete dependency provenance, tamper resistance, or cryptographic
attestation. Artifact hashes identify recorded files; they do not establish
scientific correctness.
