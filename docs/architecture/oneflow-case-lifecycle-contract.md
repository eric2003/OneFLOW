# OneFLOW Case Lifecycle and Incremental Verification Contract

## Purpose

A OneFLOW Case is a versioned, incrementally deliverable workflow, not a directory
that is considered valid only after every possible downstream operation exists.
The same lifecycle model applies to structured and unstructured meshes, imported
STL/STEP geometry, generated geometries, mesh algorithms, Python references,
native C++ implementations, format conversion, and CFD consumers.

The Case Catalog remains the discovery index. A Case may optionally reference a
separate lifecycle contract. Existing catalog entries without that field retain
their current behavior during migration.

## The lifecycle is a graph, not a fixed completion checklist

The common stage vocabulary is:

1. `geometry` — define, import, or construct source geometry / mesh-layout input.
2. `surface_discretization` — create or import boundary entities and surface mesh.
3. `volume_discretization` — create structured or unstructured volume mesh.
4. `mesh_quality` — validate topology, orientation, quality, and required regions.
5. `format_conversion` — convert between supported mesh representations (for example MSH → CGNS).
6. `oneflow_ingestion` — read the converted mesh into OneFLOW's native grid model.
7. `cfd_smoke` — execute a bounded solver smoke test and verify runtime artifacts.
8. `cfd_validation` — perform physics/numerics validation against a reference, theory, or experiment.

These are reusable stage types, not mandatory stages for every Case. A structured
mesh Case may not need geometry import or STL surface reconstruction. A mesh
reader test may stop after ingestion. A geometry-only learning Case may stop
after surface discretization. A Case declares the stages relevant to its current
scope and the dependencies between them.

## Two independent status dimensions

Every declared stage has an implementation status and a verification status.
Do not collapse them into a single green/red state.

Implementation status:

- `planned`: contract is known; implementation has not started.
- `in_progress`: implementation is being developed.
- `implemented`: implementation exists, but acceptance has not yet been established.
- `not_applicable`: this stage is intentionally outside this Case's scope.
- `blocked`: implementation is blocked by a recorded dependency or issue.

Verification status:

- `not_run`: no applicable verification has run.
- `running`: verification is currently executing.
- `passed`: all acceptance checks declared for this stage passed.
- `failed`: at least one declared acceptance check failed.
- `not_applicable`: the stage is not in this Case's scope.

A stage must not be reported as verified merely because its implementation exists.
Likewise, a later stage being planned must not invalidate tests for an earlier
implemented stage. The current deliverable is the set of stages explicitly
selected by the contract's `verification_scope`, together with their dependency
closure. Planned downstream stages are visible as planned, not treated as failed.

## Contract shape

A Case lifecycle contract uses schema `oneflow.case-lifecycle`, version 1:

- `case_id`: stable ID matching the Case Catalog.
- `contract_version`: revision of this Case's own contract.
- `verification_scope`: named milestone currently being delivered.
- `stages`: stage records with stable IDs, stage type, dependencies, status,
  implementation variants, inputs, outputs, and acceptance checks.
- `implementations`: zero or more language/tool variants (for example
  `cpp`, `python`, `gmsh`, `oneflow`). A Case need not have both C++ and
  Python on day one; absence is explicit and does not imply success.
- `acceptance`: stable test/verification IDs and the observable expected result.

The stage graph must be acyclic. Dependencies must reference declared stages.
A stage marked `passed` must be implemented and every dependency in its
dependency closure must itself be passed or explicitly not applicable.
A not-yet-implemented downstream stage may depend on a passed upstream stage;
this is the normal incremental-development state.

Inputs and outputs are stage-local contracts. Paths are relative to the stage
workspace or are logical artifact names; they must not embed machine-specific
absolute paths. Output records should identify format, dimension/topology where
relevant, and the producer stage. Boundary/region semantics must be represented
as named entities rather than inferred from filenames or hard-coded inlet /
outlet counts.

## Incremental acceptance policy

Each change is tested against the contract revision it changes:

- Adding a geometry implementation: run geometry contract tests (determinism,
  units/parameters, topology, provenance, region identity as applicable).
- Adding Python or C++ parity: run each implementation's tests and compare
  common semantic outputs, not necessarily byte-identical files.
- Adding a surface or volume mesher: test its artifact contract, mesh format,
  element families, required physical groups, and quality limits.
- Adding a converter: test the converter stage using stable fixtures and the
  current upstream artifact contract.
- Adding OneFLOW ingestion: test native grid counts, dimensions, boundaries,
  and runtime report.
- Adding CFD smoke: test process completion, bounded step count, finite
  residuals, and required result artifacts.
- Adding physical validation: add independent numerical/physical references
  and tolerances; a smoke test is not a physical validation.

The CI gate for a Case runs all currently required acceptance checks in its
declared verification scope. It does not require future planned stages. Existing
tests must remain green when a new stage is added; the new stage adds its own
tests rather than weakening or replacing upstream gates.

## Cross-language and workflow reuse

Python and C++ are implementation variants of the same stage contract. They
must share logical inputs, parameter meaning, output artifact semantics, region
names, and validation expectations. Their internal algorithms may differ.
Common converters, OneFLOW ingestion, and CFD consumers are workflow services
and should be referenced by Cases rather than copied into per-geometry scripts.

A Case's GUI entries are projections of this contract:

- Development Cases expose the available native C++ and Python implementations
  separately.
- Mesh-to-CFD Validation offers only the workflow stages currently executable
  for that Case and clearly indicates unavailable downstream stages.
- Reports preserve the full lifecycle snapshot, including every stage's
  implementation and verification status, implementation variants, current
  verification scope, and artifact hashes / lineage where available. A
  milestone report must not omit planned downstream stages.
- A OneFLOW executable is required only when the selected workflow includes
  OneFLOW conversion, ingestion, or solver stages; producer-only Python/Gmsh
  milestones remain independently runnable.
- The GUI must not infer completion from a catalog label such as
  `mesh_validation: cfd_smoke`; the lifecycle contract and actual run evidence
  are authoritative.

## Migration and compatibility

Adopt this contract incrementally. Do not migrate every historical Case at once.
Cases without a lifecycle reference continue to use the existing catalog rules.
New multi-stage Cases and Cases undergoing substantial changes should add a
lifecycle contract first. Once the loader, GUI, and CI consume the contract,
existing entries can be migrated family by family.

The 90-degree pipe bend is an initial migration target, not a special framework.
Its Python geometry/mesh path, future native C++ path, conversion path, OneFLOW
ingestion, and CFD checks are separate stages. S-shaped pipes, branched ducts,
STL imports, structured grids, and unrelated mesh algorithms use the same
lifecycle vocabulary with different stage graphs and artifacts.
