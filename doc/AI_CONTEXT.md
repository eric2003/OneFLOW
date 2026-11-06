## 2026-11-06 continuation: Case contract execution evidence linkage

The user reported 1259/1259 tests passing before this slice. That is a
pre-change baseline only; do not claim the new slice passed locally.

The verified `hexin` HEAD at the start of this slice was
`3809c0bce0ffc0bd0e0c934ea33e4db502ec7294`.

Added `capture_case_execution` as the Case-domain counterpart of
`capture_workflow_execution`. It accepts an explicit residual-runner style
run payload (command, return code, files_ok / residual_ok /
required_outputs_ok, working directory, timestamps, optional artifacts and
residual baseline path) and appends a single schema-v1 ExecutionRecord.

Contract links are accepted only when the caller passes `contract_ids`
explicitly. Matching case names, directory paths, or residual baseline
entries never create a contract link. A process exit code of zero does not
override a failed file or residual comparison.

Focused unit coverage lives in
`tools/oneflow_gui/test_case_execution_capture.py` and exercises:

- successful launch with explicit contract ID and hashed artifacts;
- refusal to record unlaunched / skipped cases;
- comparison failure recorded as failed even when exit code is zero;
- empty contract_ids remaining empty (no name/path inference).

This slice does not yet wire a Contract Workbench "Run in case runner..."
action or modify `test/test.py`. The capture function is the durable
evidence boundary; UI and residual-runner integration follow after local
validation of this capture contract.

Changed files:

- `tools/oneflow_gui/execution_capture.py`
- `tools/oneflow_gui/test_case_execution_capture.py` (new)
- `tools/oneflow_gui/test_runner.py` (register Case execution evidence tests)
- `doc/AI_CONTEXT.md`

Focused validation required:

```
python -m unittest tools.oneflow_gui.test_case_execution_capture
python -m unittest tools.oneflow_gui.test_test_inventory
python -m unittest tools.oneflow_gui.test_workflow_execution_capture
python -m unittest tools.verification.test_contract_catalog
```

No local full suite was run by the assistant. The 1259/1259 result
predates these changes. The slice remains unverified until focused tests
and full Windows/Debug CTest pass.

Next action after validation: wire an explicit Case Contract Workbench
action that launches the residual runner (or reuses the existing case
database runner) and passes the stable Case contract ID into
`capture_case_execution`. Do not infer the link from case names or paths.

---

## 2026-11-02 continuation: workflow execution evidence linkage

The user reported 1156/1156 tests passing before this slice. That is a
pre-change baseline only; do not claim the new slice passed locally.

The verified `hexin` HEAD at the start of this slice was
`d48264bfb53d2b0d442529f58d4c2be220bbdb9a`.

Added an explicit `Run in workflow runner...` action to the Workflow
Contract Workbench. It loads the checked-in RAE2822 workflow definition and
passes the stable contract ID to the existing WorkflowDialog. It does not
auto-run: the user must choose the executable/base directory and press Run.

The workflow process worker now retains per-step command, source directory,
working directory, timestamps, exit code, and declared output paths.
`capture_workflow_execution` appends immutable schema-v1 ExecutionRecords
for launched processes. It links the contract only when the runner was opened
from that contract action. Blocked or pre-launch steps remain in
WorkflowRunRecord history and are not fabricated as process runs. A process
exit code of zero does not override a failed declared-output check.

The existing Contract Workbench Verification view already separates linked
from unlinked execution history and displays revision freshness. Each step
record now carries the parent workflow run ID as an additive optional
invocation_id, preserving stable test_id semantics and backward compatibility.
The local record is not a signed SLSA provenance attestation. This slice
connects actual Workflow runs to that existing view; it does not add a second
workflow runner.

Changed files:
- `tools/oneflow_gui/execution_capture.py`
- `tools/oneflow_gui/workflow_dialog.py`
- `tools/oneflow_gui/contract_workbench_dialog.py`
- `tools/oneflow_gui/test_workflow_execution_capture.py` (new)
- `tools/oneflow_gui/test_contract_workbench_dialog.py`
- `tools/oneflow_gui/test_test_inventory.py`
- `tools/oneflow_gui/test_runner.py` (explicit Python unittest registry)
- `.github/workflows/python-tools.yml`
- `doc/plans/verification-contract-workbench.md`
- `doc/source/workbench_architecture.rst`
- `tools/verification/README.md`

Focused validation required:
`python -m unittest tools.oneflow_gui.test_workflow_execution_capture`
`python -m unittest tools.oneflow_gui.test_contract_workbench_dialog`
`python -m unittest tools.oneflow_gui.test_workflow_repository`
`python -m unittest tools.oneflow_gui.test_workflow_graph_runner`
`python -m unittest tools.oneflow_gui.test_test_inventory`

No local test suite was run by the assistant. The 1156/1156 result
predates these changes. The slice remains unverified until Windows focused
and full tests pass.

Next action after validation: add comparable explicit execution linkage for
Case contracts, using the existing case database/regression runner records.
Do not infer a contract link from matching case names or paths alone.

---

## 2026-11-02 continuation: first typed Workflow contract

The user reported 1156/1156 tests passing before this slice. This is the
pre-change baseline only; no local test result is claimed for the changes
below.

Inspected the actual `hexin` HEAD:
`e8369e5b09cfab84cb88c324df3b62a6b33731df`.

Added the first Workflow-domain contract using the checked-in
`tools/oneflow_gui/examples/rae2822_multicase.json` example. The typed payload
records the mesh producer, Roe-SA consumer, direct dependency, stable artifact
ID, producer output path, and consumer `gridFileName` parameter binding.

The Qt-independent catalog now validates safe relative paths, unique step and
artifact IDs, producer resolution, direct dependency edges, and consistency
between the typed payload and the checked-in workflow definition. It does not
duplicate topological sorting or runner behavior. Those remain owned by the
existing workflow planning and graph execution services.

The Contract Workbench domain hub now exposes Workflow. Its Workflow Graph tab
shows the declared steps and artifact bindings and provides read-only access to
the source definition. Linked Python plan/runner tests remain
`linked_not_executed`; no workflow or CFD execution is inferred.

Changed files:
- `tools/verification/contract_catalog.py`
- `tools/verification/test_contract_catalog.py`
- `tools/verification/contracts/workflow/rae2822_grid_to_solver.json`
- `tools/oneflow_gui/contract_workbench_dialog.py`
- `tools/oneflow_gui/test_contract_workbench_dialog.py`
- `tools/oneflow_gui/test_test_inventory.py`
- `doc/plans/verification-contract-workbench.md`
- `doc/source/workbench_architecture.rst`

Focused validation required:
`python -m unittest tools.verification.test_contract_catalog`
`python -m unittest tools.oneflow_gui.test_contract_workbench_dialog`
`python -m unittest tools.oneflow_gui.test_workflow_plan`
`python -m unittest tools.oneflow_gui.test_workflow_graph_runner`
`python -m unittest tools.oneflow_gui.test_test_inventory`

The assistant has not run the Windows test suite. The 1156/1156 result
predates this slice. The workflow contract is implemented but unverified until
these focused tests and the full suite pass on the resulting revision.

Next single action: validate this slice on Windows, then add execution-record
linkage only when a real workflow run can retain exact revision, command,
configuration, outcome, and artifact references. Do not add a generic runner.

---

# AI and Developer Project Context

## 2026-11-01 continuation: first typed CFD Case contract

The user reported 1126/1126 tests passing before this slice. That is the
baseline only; no local test result is claimed for the changes below.

The active architecture workstream remains the cross-domain Contract Workbench.
The repository-root `gui/` directory remains temporarily out of scope. The
latest hexin HEAD before handoff is
`e0d315a46f63c49a7af0487a753b95a63739219b`.

Added `tools/verification/contracts/case/laminar_flat_plate_slau2.json`,
a real Case-domain definition for `test/cfd/plateuns2dslau2`. It is already
selected by `test/suites/cpu-serial.txt`. The definition records:
- 2-D laminar compressible Navier-Stokes, `NsSolver`, SLAU2, 50 steps;
- control scripts and grid inputs;
- the exact five paths from `autotest/compare.txt`;
- `results/res.dat` as the residual output;
- `test/baselines/residual-baseline.json` with absolute tolerance 1e-15;
- retained Blasius Cf/u/v references, explicitly marked as not automatically
  compared by the current runner.

`contract_catalog.py` now has immutable `CaseContractData` as a typed
domain-specific payload attached optionally to `ContractSummary`. It checks
that case paths and input/reference files exist, the case is selected by the
declared suite, the output list exactly matches `autotest/compare.txt`, and
the residual database:
- has the expected schema;
- contains the named case;
- contains every declared residual output;
- contains iteration columns `iter` and `sub-iter`;
- has the configured number of residual rows;
- uses the declared absolute tolerance.

Tests cover Case payload adaptation, missing payload, combined catalog
discovery, tolerance drift, and residual-history length drift. The catalog
does not launch CFD and does not mark the case as executed. A baseline match
is regression reproducibility, not independent physical validation. Blasius
reference files remain declared-only evidence until an explicit comparison
is implemented and scientifically reviewed.

Documentation updated:
- `tools/verification/README.md`
- `test/README.md`
- `doc/plans/verification-contract-workbench.md`

Important commits in this slice:
- typed payload: `253c51948cf94e88f96ca71deb68795a1574e043`
- Case definition: `c6410ec05bb4558b7db637faa43d7e27a22f4450`
- catalog checks: `80f4c46aa62237c5b45628ed384095b074e22a5d`,
  `3c8882f4e3afd2f685f2ef06dab539bd98504098`,
  `811b36cbb6451bcaed0c71debc1152d2aeb64320`
- test coverage: `5b6d854085e0efbefbe71f8527518300a0b830ad`,
  `6607bb791706765cf81b6da39ff3ef76aab2efc4`,
  `177e7e2b9688f29427ee4a808bbefc6736b2ebe7`,
  `be6ea51bb09e0fc4cd91e1401e07d0394de88ae1`
- documentation: `bb4300822edae9f3286aae5b906f176347947327`,
  `2c3125f4c943de8c44af672d242cc59b3a4cd4e8`,
  `a7c5451343d924c93605efe848443b2925004f47`,
  `003947f073853491983752c8f4bd69de2f1d9fc4`,
  `e0d315a46f63c49a7af0487a753b95a63739219b`.

Focused validation required:
`python -m unittest tools.verification.test_contract_catalog`.
The previous 1126/1126 result predates this slice. Next, run focused tests
and full Python tests. Once validated, add a second Case/CFD contract only if
it exercises a genuinely different payload (e.g. turbulence or grid-only
conversion), then use a Workflow contract to test the same shared evidence
model. Do not add a generic runner or GUI.

## 2026-11-01 continuation: Contract Workbench becomes cross-domain data

The user reports 1126/1126 tests passing before this slice. That result is the
baseline only; it does not validate the changes below.

The Contract Workbench direction is an active architecture workstream, not a
GUI-only feature. The repository root `gui/` remains out of active scope, and
do not spend effort expanding `tools/oneflow_gui` while the C++ and contract
foundations need attention.

Implemented the first declared Runtime contract:
- `tools/verification/contracts/runtime/process_lifecycle.json`
- `tools/verification/contract_catalog.py`: `ContractSummary` now exposes
  common architecture/status fields; `adapt_declared_contract` loads schema
  version 1; `load_contract_catalog` discovers declared contracts alongside
  Numerical and Geometry compatibility adapters.
- `tools/verification/test_contract_catalog.py`: tests the declared runtime
  contract, schema-version rejection, evidence status, and reference validity.
- `tools/verification/README.md` and
  `doc/plans/verification-contract-workbench.md` describe the vertical slice.

The runtime contract links exact production implementation paths and three
GoogleTest declarations. They remain `linked_not_executed`; catalog loading
checks source declarations but never infers a pass. Contract status
`implemented` is not a verification result. Its claims explicitly exclude
physical CUDA/HIP driver behavior, multi-rank MPI teardown, all solver resource
release paths, and CFD numerical correctness.

Commits:
- catalog model: `ffe646c4bfb82d33f3b5121d8eeb0d3c2d9d4aa7`
- runtime contract definition: `76934c00585ed7245035ffe674505afcf39f4546`
- tests: `170293fbd9f90269df2cf3882efb435a09be119d`
- verification README: `3ed2596e16d409a7db7dc23458d1770fd79c70eb`
- architecture proposal: `7fec4524fd6df556de8f445543becebbc94970be`

Focused local validation required:
`python -m unittest tools.verification.test_contract_catalog`.
No assistant-side Windows build or test result is claimed.

Next architectural decision: use a second genuinely different domain
(prefer a Case/CFD regression contract, with workflow as another candidate) to
test whether the shared envelope is sufficient. Do not introduce a generic
runner, plugin framework, or wholesale JSON migration yet. Case definitions
should normally remain evidence fixtures linked to multiple contracts, not
become one contract per case. Keep claims, evidence scope, exact command,
revision, backend/platform, outcome, and artifacts distinct. Once two real
domains exercise the same extension point, consider a typed evidence payload
or adapter interface.

## 2026-11-01 continuation: production process-level accelerator retry

The user reports all 41 registered Python command groups passed. The root
`gui/` directory is temporarily out of active development scope; do not spend
the next slices expanding its legacy test inventory. Keep the existing
`tools/oneflow_gui` infrastructure stable, but prioritize C++ core,
numerical algorithms, runtime architecture, documentation, and CFD regression.

Added a production-linked lifecycle test:
`SimuImpCaseLifecycleIntegration.FailedAcceleratorFinalizationKeepsProcessRuntimeActiveForRetry`.
It injects a retryable backend under the HIP backend kind, forces the first
`AccelRuntime::Finalize()` call to throw through the real
`SimuImp::FinalizeEnvironment()` entry point, and verifies the process
lifecycle remains Active and the backend remains initialized. A second
finalization succeeds and commits the lifecycle transition.

Design record: `doc/plans/accelerator-runtime-shutdown-order.md`, section 8.
Implementation commit: `5e38b60444c6090e372a77326e987bfdad4da99d`.
Documentation commit: `f111469cac3458701f9a1a2a2e9ef4902db6d75e`.

Status: implemented, not yet validated on the user's Windows build. The fake
HIP-kind backend verifies process-level retry ordering without requiring HIP
hardware; it does not replace HIP driver or multi-rank MPI validation.

Next validation: build and run the focused
`SimuImpCaseLifecycleIntegration.FailedAcceleratorFinalizationKeepsProcessRuntimeActiveForRetry`,
then `AccelRuntimeTest.FailedFinalizationPreservesBackendForRetry`, followed
by full Debug CTest. If they pass, proceed to the remaining case-teardown
failure boundaries and representative CFD regressions.

Last updated: 2026-11-01
Working branch: `hexin` only. Latest inspected baseline for the shutdown-order slice: `6e550260578a7ad3490793791f9a26feb1e670ae`. Current hexin HEAD after fixture correction and handoff updates: see the latest commits in this section.
User-reported local validation on the pre-change baseline: **1123/1123 tests passed**. This does not validate the shutdown-order change added afterward.

## 2026-11-01 continuation: lifecycle integration test fixture correction

The user reported 1123/1124 tests passing after the accelerator shutdown-order
change. The sole failure was
`SimuImpCaseLifecycleIntegration.ProcessFinalizationReleasesContextBackendStatesBeforeRuntime`,
which threw in `Prj::ResolveCaseDir()` during `SimuImp` construction. The
test passed the relative placeholder `unused-case` before initializing
`Prj::current_dir`; therefore it never reached process runtime setup or the
shutdown-order assertions.

Corrected the test fixture on `hexin`:
- Call production `Prj::Init()` before constructing `SimuImp`.
- Reuse the existing RAII `TemporaryCase` fixture and pass its absolute
  directory to the constructor.
- Document why the test requires an initialized project environment and an
  absolute case path.

Updated design record:
[Accelerator runtime shutdown order](plans/accelerator-runtime-shutdown-order.md).

Correction commits:
- Test fixture: [aeae2fb](https://github.com/eric2003/OneFLOW/commit/aeae2fb51e719d0d6a1a8a9786f22f18dfb9230b).
- Design/validation record: [298e022](https://github.com/eric2003/OneFLOW/commit/298e022c8feddf4fea1f46fec05ccca4031e74d5).

Status: **Fixture corrected; rerun not yet verified**. The 1123/1124 result
belongs to the pre-correction revision. First rebuild and run the focused
integration test, then the full Debug CTest suite. If the focused test passes,
inspect any remaining failures before proceeding to the next runtime
finalization exception-safety slice.

## 2026-11-01 continuation: retryable accelerator shutdown

The user reports **1124/1124 tests passed** on the preceding hexin revision
after the process-finalization test fixture was corrected. This is the
pre-change baseline for the current slice.

Starting revision: `ba067226fcfefad896648d24bd7daa75b49b3a20`.

Source review found a dependency bug in the previous shutdown sequence:
`ExecuteAllCleanupSteps(ClearAccelStates, ReleaseShared, FinalizeAccelRuntime,
HXFinalize)` treats all four operations as independent. They are not.
`HXFinalize()` must not run when accelerator finalization fails, because the
backend may still require the process/parallel runtime. The previous helper
would continue into HXFinalize after a failed FinalizeAccelRuntime.

Implemented on `hexin`:
- Split independent pre-runtime releases from dependent runtime finalization
  in `SimuContext::FinalizeEnvironment()`.
- Clear context-owned accelerator states and release the shared HIP backend
  together; attempt both, then stop if either reports failure.
- Call `FinalizeAccelRuntime()` separately. If it throws, skip
  `HXFinalize()` and leave the process lifecycle phase Active.
- Call `HXFinalize()` and mark the process finalized only after accelerator
  shutdown succeeds.
- Documented the backend contract: if `AccelBackend::Finalize()` throws,
  the backend must remain valid and safe to retry. `AccelRuntime::Finalize()`
  already retains its backend unique_ptr and initialized flag when the virtual
  call throws.
- Added `AccelRuntimeTest.FailedFinalizationPreservesBackendForRetry`, which
  injects one backend finalization failure, verifies the backend remains
  accessible/initialized, then verifies a successful retry.

Design and architecture records:
- [Accelerator runtime shutdown order](plans/accelerator-runtime-shutdown-order.md)
- [Runtime architecture](source/runtime_architecture.rst)

Industry basis: C++ Core Guidelines E.4/E.17 (invariants and exception
propagation) and C.36/C.37 (non-throwing destructors). The important design
distinction is between independent cleanup steps, which should all be
attempted, and dependent lifecycle transitions, which must stop when a
precondition has not been established.

Status: **Implemented, unverified**. The user's 1124/1124 result predates
these changes. Run the focused accelerator runtime test, the production
SimuImp lifecycle integration test, and full Debug CTest on the exact updated
hexin HEAD. Also check Windows/Linux and documentation CI. HIP hardware failure
behavior and multi-rank MPI remain separate validation items.

## Current follow-on: accelerator runtime shutdown ordering

The user reports **1123/1123 tests passed** on baseline hexin HEAD
6e550260578a7ad3490793791f9a26feb1e670ae. The current slice addresses a lifecycle
gap discovered while tracing SimuContext, EulerDomainStateRegistry, DeviceBuffer,
and AccelRuntime ownership.

- SimuContext owns EulerDomainStateRegistry; registered backend states may own
  device allocations.
- TeardownCase() already clears context states and releases the shared HIP flux
  backend before process finalization.
- FinalizeEnvironment() was independently callable but did not clear any
  remaining context-owned states before FinalizeAccelRuntime(). A caller that
  skipped TeardownCase() could therefore destroy backend state after runtime
  shutdown.
- FinalizeEnvironment() now clears the registry and releases HipFluxBackend's
  shared instance (when HIP is enabled) before FinalizeAccelRuntime(), then
  calls HXFinalize(). This is defensive process cleanup, not a replacement for
  full TeardownCase() of solver/global case state.
- A production-linked integration test inserts a state whose destructor
  records whether AccelRuntime is still initialized, then finalizes without
  TeardownCase() and checks destruction order.
- Design record: [Accelerator runtime shutdown order](plans/accelerator-runtime-shutdown-order.md).
- Runtime architecture guide updated in source/runtime_architecture.rst.

Industry references: C++ Core Guidelines R.1, C.31, C.36, C.37, E.4 and E.7.
The important distinction is that RAII manages each object's lifetime, while
the owning context must still order that lifetime before the allocator runtime
it depends on.

Status: **Implemented, unverified**. The user-reported 1123/1123 result predates
this slice. Build/run the production integration target, full Debug CTest,
and representative CPU regressions on the exact latest hexin HEAD. Also check
Windows/Linux CI. HIP-enabled compile/runtime coverage remains distinct from
the CPU ownership-ordering test.


## Previous mainline: R1 case configuration rollback and lifecycle

The mainline is to make the production case bootstrap exception-safe without confusing configuration rollback with full solver teardown.

Production path: `SimuImp::Run() -> RunCase() -> PreProcess() -> InitSimu() -> SimuContext::SetupCaseEnvironment() -> InitializeCaseConfiguration(caseDir) -> ReadControlInfo(caseDir) -> ReadCaseIOSettings() -> ApplyCaseIOSettings()`.

Current implementation:
- `codes/main/src/SimuContextCaseEnvironment.cpp` now wraps the production configuration bootstrap in `SimuCaseLifecycleDetail::ExecuteWithRollback()`.
- On bootstrap failure, `SimuContextDetail::ClearCaseConfigurationState()` clears partial DataField/DataPara state, the case log binding, and the project-root binding. The rollback helper preserves the original bootstrap exception if rollback itself throws.
- `IsEnvironmentReady()` remains false until the entire bootstrap succeeds.
- `SimuImp::RunCase()` still owns full `PostProcess() -> TeardownCase()` cleanup for solver/runtime resources on both setup and task failures. The new rollback is a narrower guarantee for callers of the public Context setup API; it does not replace full teardown.
- `tests/io/control_configuration_ingestion_test.cpp` now exercises the production Context setup path: the first script commits a value, the second script is missing, and the test verifies the partial value and project binding are cleared. It then creates the missing script and retries successfully on the same context.
- The separate two-case bootstrap test covers configured roots on one case and default-root rebinding on the next. It still uses the production configuration composition and real control reader.

System-level reason: case configuration now has a strong failure boundary independent of whether the caller is `SimuImp` or another Context consumer. This prevents partially ingested global configuration from leaking after a failed bootstrap. It does not yet prove full solver/device resource release through the actual `SimuImp::RunCase()` callback wiring.

## Validation status

- User reports 1096/1096 local tests passed on the pre-change baseline. No local build or CTest run is available in this environment for this slice.
- Runtime implementation commit: [6a325f7](https://github.com/eric2003/OneFLOW/commit/6a325f71ede48251df3d50e9e65d3a98f8b5dafe).
- Integration assertions and retry test: [4de3f86](https://github.com/eric2003/OneFLOW/commit/4de3f869affe58a5ac495a4cb90f12024bb3a6f7), [729811b](https://github.com/eric2003/OneFLOW/commit/729811b0996c929ec9b566bee643ab13b0897b75).
- Runtime architecture guide update: [abf2bd2](https://github.com/eric2003/OneFLOW/commit/abf2bd241d2139177630a2a99e7bb21d531e591a).
- Exact-HEAD Sphinx check succeeded; Windows and Linux CI were in progress at the last poll: [Sphinx](https://github.com/eric2003/OneFLOW/actions/runs/41280382266), [Windows](https://github.com/eric2003/OneFLOW/actions/runs/41280382252), [Linux](https://github.com/eric2003/OneFLOW/actions/runs/41280382250).
- This environment has no mounted OneFLOW Windows checkout, so no local VS2022/MSBuild/CTest result can be produced here.
- C++ source files remain ASCII-only by default; Python and mathematical presentation may retain useful Unicode characters.

## Remaining R1 gap

R1 remains **Partial**. The production integration target exercises the real `SimuImp::RunCase() -> PostProcess() -> TeardownCase()` wiring after a setup failure. A new success-path integration test is now implemented but unverified: it runs two separate `SimuImp` cases through the production `RunCasesWithProcessLifecycle()` orchestration, overrides `ToyModel` with a marker task, and checks case-state cleanup plus process-runtime continuity between cases. Remaining gaps are injected failures at individual solver/device release points, multi-rank MPI behavior, and representative full CFD regressions. The shared `OneFLOW_RuntimeObjects` CMake object collection now lets the executable and integration target consume the same production translation units without creating a static-library boundary. Do not substitute fake teardown and call it end-to-end.

The earlier investigation record remains `doc/plans/setup-case-environment-integration-test.md`. Its original findings are historical; this handoff and `doc/source/runtime_architecture.rst` record the newer rollback behavior.

Next single action: build and run the focused production lifecycle integration target on the exact latest commit, then run the full Debug CTest suite and representative grid/CFD regressions. If those pass, update the R1 evidence; keep teardown fault injection and multi-rank MPI as explicit remaining gaps. Detailed findings are in `doc/plans/setup-case-environment-integration-test.md`.
## Previous follow-on: accelerator buffer ownership

The latest branch focus has moved to accelerator resource reliability after the recent
transactional backend initialization work. Inspected starting HEAD: d9e25a7b671e6104c0d116ae5c561e5d12ba09d0 (Register accelerator runtime tests). The user reports 1119/1119 local tests passed on that pre-change baseline.

- DeviceBuffer<T>::Resize() previously leaked its replacement allocation if device-to-device copy failed.
- Its destructor and noexcept move assignment called a backend deallocator that may throw.
- The implementation now guards candidate allocations, preserves the old buffer on allocation/copy failure, rejects byte-count overflow and out-of-capacity host copies, commits replacement ownership before releasing the old allocation, and uses non-throwing cleanup from destruction/move assignment.
- Explicit Reset() reports release errors but clears object ownership first. Buffers must still be released before accelerator runtime finalization; a destructor cannot free through a backend that has already been finalized.
- Focused fake-backend tests cover failed growth cleanup/original-buffer preservation and non-throwing destruction.
- Design record: [DeviceBuffer resource lifecycle](plans/device-buffer-resource-lifecycle.md).

Implementation: [9728e46](https://github.com/eric2003/OneFLOW/commit/9728e46f97204a04a0bc188e770fde2e3d453e23). Focused test coverage and tracking backend: [df8dd00](https://github.com/eric2003/OneFLOW/commit/df8dd0008356e85ff0021cbb221535d9045069fd), [3bd2bf8](https://github.com/eric2003/OneFLOW/commit/3bd2bf84ca4f7b2abd71be94e2e77faf513ffffa), [e6d3e90](https://github.com/eric2003/OneFLOW/commit/e6d3e902e2f64dacb9ec1cb9a27cfbc3ca4cd2e5), [c19ef29](https://github.com/eric2003/OneFLOW/commit/c19ef29558f5a62ebe34e595d40976ad2af59aac), [a48c2b6](https://github.com/eric2003/OneFLOW/commit/a48c2b6f79651ec24e98d979a6a900d61a43dbf8), [74de922](https://github.com/eric2003/OneFLOW/commit/74de9222c90c7bce19602c209d4ddbdd1ab10b21).
Status at the time of that slice: implemented. The subsequent 1123/1123 local test report covers the DeviceBuffer changes, not the newer shutdown-order change. No local compiler/CTest execution is available in this environment. Next action is to build and run accel_runtime_test, then the full Debug CTest suite and CPU regression cases on the exact updated HEAD.

## 1. Read this first

For a quick restart, follow this order:

1. **This page** — repository constraints, current focus, known gaps, and the next useful action.
2. [Whole-system architecture review](source/architecture_review.rst) — common CFD/CAE concerns, OneFLOW's cross-system risks, and the required zoom-out review at milestones.
3. [Runtime architecture and development guide](source/runtime_architecture.rst) — source-backed entry points, runtime lifecycle, configuration, and current boundaries.
4. [Runtime refactoring roadmap](plans/runtime-refactoring-roadmap.md) — staged plan, acceptance criteria, and per-slice findings.
5. [Documentation consolidation plan](plans/oneflow-documentation-consolidation-plan.md) — canonical documentation locations and maintenance rules.
6. The relevant focused tests and source files named below; then check CI for the exact commit being assessed.

Do not start by reading the whole repository. Use the architecture map to identify a concrete path, then inspect its callers, implementation, tests, and evidence.

## 2. Repository and collaboration constraints

- Repository: [`eric2003/OneFLOW`](https://github.com/eric2003/OneFLOW).
- Default working branch for this continuing effort: `hexin`. Base changes on the latest `hexin` HEAD, not `master`, unless the user explicitly changes this instruction.
- Keep changes merge-compatible with `hexin`; do not push or submit this work to `master`.
- Prefer small, evidence-backed changes that solve a real behavior or maintenance problem. Do not refactor merely to introduce a pattern.
- Follow direct call paths and state/resource ownership before changing interfaces.
- Preserve compatibility unless an intentional behavior change is stated and tested.
- Keep C++ source files ASCII-only.
- Do not claim a test ran, passed, or a CI job succeeded without evidence for the exact relevant revision.
- Keep public documentation aligned with current behavior. Distinguish implemented behavior, partial work, planned work, and unverified assumptions.
- The user values architecture discussion, meaningful tests, concise reasoning, and exact commit links.

## 3. Whole-system review checkpoint

The published [architecture review](source/architecture_review.rst) is now the common checklist for substantial AI-assisted changes. It deliberately separates general CFD/CAE concerns from verified OneFLOW behavior. Keep these cross-system risks visible while working locally: implicit global runtime state, configuration parsing versus validation/application, mutable mesh data and derived connectivity, workflow artifact provenance and execution locations, repeated/failed-run cleanup, and the gap between documented plans and actual test/CI evidence. At each meaningful milestone, state which system-level concern improved, what remains untouched, and whether the priority order should change.

## 4. Current focus: R1 case configuration and lifecycle

Source and tests:
- `codes/main/src/SimuContextCaseEnvironment.cpp`
- `codes/main/src/SimuContextBootstrap.cpp`
- `codes/main/src/SimuContextCaseState.cpp`
- `codes/main/include/SimuCaseLifecycle.h`
- `tests/io/control_configuration_ingestion_test.cpp`
- `tests/main/simu_case_lifecycle_test.cpp`
- `tests/main/simu_imp_case_lifecycle_integration_test.cpp`
- `doc/source/runtime_architecture.rst`
- `doc/plans/runtime-refactoring-roadmap.md`

Verified contracts:
- Case setup publishes readiness only after the complete production bootstrap succeeds.
- A bootstrap failure rolls back partial DataField/DataPara state, case log binding, and project-root binding; the original setup exception is preserved.
- The same Context can retry after the failed manifest is repaired.
- `SimuImp::RunCase()` remains responsible for full case teardown; process finalization remains a separate owner.

Verified by production-linked integration: a real setup failure through `SimuImp::RunCase()` resets context-owned case state and leaves process runtime active until explicit finalization. Implemented, unverified: a successful two-case sequence using separate production `SimuImp` instances, checking per-case reset and process runtime continuity. Still open: injected teardown failures in the concrete release chain, multi-rank MPI, and full numerical regressions.
## 5. Next actions, in priority order

1. Check exact-HEAD Windows/Linux CI for the current `hexin` revision; distinguish queued, in-progress, and completed jobs.
2. Build and run the new successful sequential-case integration contract on the exact commit.
3. Run the full Debug CTest suite and representative grid/CFD regression cases.
4. If successful, update R1 evidence; then investigate teardown fault injection at concrete solver/device release points.
5. Keep multi-rank MPI validation as a separate acceptance item and reassess R1 status only against the required evidence.

Suggested next mainline: validate production sequential-case isolation, then close the remaining concrete cleanup-failure and MPI gaps—not another generic lifecycle helper or unrelated local refactor. The whole-project investment ranking is documented in `doc/source/architecture_priority_roadmap.rst`.
## 6. Status vocabulary

Use these labels consistently in this page and related records:

- **Verified** — confirmed by source inspection, a named test, or a successful CI check; state which evidence.
- **Implemented, unverified** — code exists, but required test/build/CI evidence is not yet available.
- **Partial** — some behavior is covered, but an acceptance criterion remains open.
- **Planned** — not implemented.
- **Blocked** — progress depends on a named missing capability, test seam, or external condition.
- **Deferred** — intentionally postponed, with the reason recorded.

A green unrelated test suite is not evidence that a changed contract is correct. A queued or in-progress workflow is not a pass.

## 7. Reusable stage records

For a substantial refactoring stage, use the [refactoring stage record template](plans/templates/refactoring-stage-record.md). Put detailed stage-specific problem statements, acceptance criteria, source paths, tests, decisions, and evidence in that record. Keep this page as a concise index of the overall current state and highest-priority next actions; link to detailed records rather than duplicating them.

## 8. How to maintain this handoff

Update this page when a milestone materially changes the state of the project, not after every command or minor edit.

For each meaningful update, record:

- inspected branch and base commit;
- the behavior or architecture question being addressed;
- current source path and relevant tests;
- what changed and what deliberately did not change;
- exact validation evidence, including failures and unavailable checks;
- unresolved gaps and the next recommended action;
- links to implementation and documentation commits.

Keep detailed reasoning and acceptance criteria in the relevant architecture page or roadmap. Keep reusable build/test commands in runbooks. Keep CI or runtime evidence in evidence records when a durable standalone record is useful. Avoid duplicating entire documents here; this page should point to them.

When starting a new AI session, first compare the current `hexin` HEAD with the implementation commit recorded above, inspect any changes since that commit, and then update this status snapshot only after verifying the new state.


## Evolutionary architecture

The new design page `doc/source/evolutionary_computation_architecture.rst` defines the long-term target: small computation units with explicit typed ports, preconditions/postconditions, side effects, lifecycle, failure policy, provenance, and verification evidence. The JSON Schema starter is `doc/specs/computation-unit.schema.json`. Both are proposals; no general plugin ABI or automatic pipeline runtime is claimed. The intended sequence is configuration ingestion first, then mesh/topology/boundary/result intermediate representations, and only later composable solver pipelines if real use cases justify them.


Latest source slice (implemented, unverified): `ReadOneFLOWScriptFile` overloads moved from `ParaFile.cpp` to `ControlScriptParser.cpp`; `GetParameterArraySize` moved to `ControlParameterResolver.cpp`. This makes the injected resolver seam in `ConfigLoader` explicit: parser tests can supply deterministic array-size resolution, while production uses the real database-backed resolver. The legacy `ParaFile.cpp` now retains inspection helpers. See `doc/plans/setup-case-environment-integration-test.md` for expected failure modes, partial-commit/transaction design question, and validation boundaries.

A focused `tests/io/control_parameter_resolver_test.cpp` and CMake target now cover literal array lengths, named lengths through the real DataBase, and missing-name failure. This is component evidence for the resolver only; full `ReadControlInfo(caseDir)` ingestion and parallel/database-dump integration remain open.


## Architecture investment and documentation policy

The ranked system-wide investment sequence is recorded in
`doc/source/architecture_priority_roadmap.rst`: (1) runtime lifecycle and
case isolation, (2) mesh core/topology/boundary lifecycle, (3) numerical
algorithm verification, (4) computation-unit contracts and dataflow
composition, (5) workflow/GUI, and (6) documentation automation. This is a
priority order, not a claim that lower-ranked work is unimportant.

Use source comments for stable API contracts and generate API reference pages
from Doxygen-compatible comments where practical. Comments should describe
intent, units, ownership/lifetime, invariants, preconditions, postconditions,
side effects, and failure semantics. Keep cross-file architecture, mathematical
derivations, numerical assumptions, experimental results, diagrams, and
reference comparisons in curated `doc/source` pages. Generated API docs
complement these guides; they do not replace them.

The architecture priority roadmap uses Sphinx math and text diagrams to make
the ranking, configuration transformation path, mesh representation layers,
and numerical verification chain explicit. Continue using equations, figures,
and diagrams where they explain actual algorithms or boundaries, and label
planned designs as proposals rather than implemented behavior.


## Documentation structure and source-comment policy

Keep `doc/plans/` separate from `doc/source/`. Plans record intended changes,
sequencing, acceptance criteria, risks, hypotheses, and unresolved questions.
The Sphinx source manual records current architecture, verified contracts,
algorithm derivations, scientific assumptions, figures, and stable design
knowledge. Promote completed plan conclusions into the relevant source page;
retain only useful decision history and remaining work in the plan. The policy
and directory map are described in `doc/README.md`.

Use Doxygen-compatible source comments for stable API contracts: purpose,
units, ownership/lifetime, invariants, preconditions, postconditions, side
effects, and failure semantics. Generate API reference pages from those
comments where practical. Keep cross-file architecture, derivations, formulas,
numerical verification, and conceptual diagrams in curated `doc/source`
pages. Source-generated API documentation complements, but does not replace,
the system and algorithm guides.

The 2026-10-30 IO adapter test adds two contracts: typed values reach the real
DataBase, and repeated declarations are committed in document order with the
later value replacing the earlier value. The adapter is incremental, not
transactional; changing that behavior requires separate compatibility analysis.

## 2026-10-30 continuation after the 976/976 baseline

The latest inspected `hexin` starting HEAD was
`d00d86c64cf940f57106a68f62f406c8b6b36fc5`. The user reports 976/976 local
tests passed before this continuation; no new local test run is available here.
Windows and Linux CI for that exact HEAD were still in progress at inspection:
[Windows](https://github.com/eric2003/OneFLOW/actions/runs/41055632631),
[Linux](https://github.com/eric2003/OneFLOW/actions/runs/41055632649). Do not
report them as passing until their exact-run conclusions are checked.

The source review reconfirmed the production reader sequence:
manifest/script parsing -> process greeting/collection -> database broadcast
(one-process broadcast is a no-op after valid parallel setup) ->
`log/database.log` dump. A real serial integration test must call
`SetUpParallelEnvironment()`; static-zero parallel fields are not a valid
fixture. The smallest checked-in manifest fixture is
`test/grid/rae2822convert` and it selects the Grid task, which must not be
executed by a setup-only test.

CMake review found that `codes/CMakeLists.txt` recursively collects the
production code into the `OneFLOW` executable with global compile definitions
and dependencies, while IO tests select individual source files. No shared
configuration-runtime library currently exists. Do not copy the whole
executable source list into a test. Map the true dependency closure first and
create a shared production target only if a cohesive reusable boundary emerges.

Two documentation commits on `hexin` record this continuation:
- `cd806307af0566a889ba89c5da2baa3f1ae7d495` corrects the adapter source
  boundary: `ControlScriptParser.cpp` invokes the adapter implemented in
  `ConfigDatabaseAdapter.cpp`.
- `baea9f4e5342c84bde7e276b21a7f9ecefaeed45` records the reader execution
  sequence, serial initialization contract, CMake target boundary, and fixture.

Next action remains a real multi-file `ReadControlInfo(caseDir)` integration
test, following dependency closure mapping. R1 remains Partial; no evidence yet
for full reader integration, effective path application, sequential-case
isolation, or setup-failure cleanup.
## 2026-10-30 continuation: manifest discovery seam

Starting from the user-reported 976/976 local baseline, the latest `hexin`
source slice is `edacedd47f38c10b114f72c330af9243e16b5b8b`. No local build or
CTest run was performed in this environment. Exact-HEAD CI is queued:
[documentation](https://github.com/eric2003/OneFLOW/actions/runs/41057437914),
[Linux](https://github.com/eric2003/OneFLOW/actions/runs/41057437918),
[Windows](https://github.com/eric2003/OneFLOW/actions/runs/41057437955).

A narrow production extraction moved both `ReadScriptFileNameList` overloads
from `ControlFileReader.cpp` into `ControlManifestReader.cpp`. The public
contract in `ParaFile.h` now states that manifest paths are appended and their
order is preserved. `tests/io/control_manifest_reader_test.cpp` exercises the
case-directory overload using a temporary manifest and verifies blank-line
handling, case-root resolution, and order; it is registered with
`TextFileParser_test`. The new translation unit is also discovered by the
recursive production source collector.

This is component-level coverage for manifest discovery only. It does not
invoke script parsing/commit, `ReadMultiScriptFiles`, parallel coordination,
`ReadControlInfo(caseDir)`, or database dumping. R1 remains Partial. The
architecture rationale and exact boundary are in
`doc/plans/setup-case-environment-integration-test.md` and
`doc/source/case_parameter_design.rst`.

Next: inspect exact-HEAD Windows/Linux CI. If the new test target builds, map
the dependency closure from `ReadMultiScriptFiles` and
`ReadControlInfo(caseDir)`; then decide whether a production-shared
configuration-runtime target is justified. Do not broaden the test target by
copying the entire executable source list.


## 2026-10-30 continuation: manifest append contract

The user reports a new local baseline of 977/977 tests passed. This is the
pre-change baseline for the current continuation; it is not evidence for the
following GitHub edits.

Starting from `hexin` HEAD
`4cb14c0b28af7dd319cdc62e23075c8df1f066d9`, the manifest reader test was
strengthened to initialize its destination vector with an existing entry,
then assert that the entry remains at index zero while the three manifest
entries are appended in source order. This pins both parts of the public
contract: append semantics and stable manifest ordering. The production
reader is unchanged. Commit:
`13ce812c50b813b94afded7dbb1257a370876eaa`.

No local compiler or CTest run was available in this environment. The test
change must be validated on the user's Windows build and exact-HEAD CI before
being treated as verified.

System-level interpretation: manifest discovery is now a better specified
boundary, but it remains only the first stage of the R1 configuration path.
Do not expand this seam into a generic configuration framework. The next
high-value step remains mapping the real dependency closure for
`ReadMultiScriptFiles` and `ReadControlInfo(caseDir)`, then exercising
production parsing, typed database commit, serial coordination, diagnostic
output, sequential-case isolation, and failure cleanup.

Documentation direction: keep source comments focused on API contracts,
units, ownership, preconditions, postconditions, side effects, and stable
algorithm identifiers. Keep derivations, equations, assumptions, diagrams,
comparisons, and curated references in `doc/source`. Doxygen API output
should complement Sphinx rather than replace curated theory and architecture
pages. Prefer a small pilot on one stable module before adding Breathe or a
larger documentation toolchain.


## 2026-10-30 continuation: standardized engineering task specification

The user reports a pre-change local baseline of 977/977 tests passed. This is
user-reported evidence for the baseline, not a post-change result.

The task-expression requirements have been consolidated into maintained
specifications:
- `doc/AI_TASK_TEMPLATE.md`: canonical English engineering task specification.
- `doc/AI_TASK_TEMPLATE_zh_CN.md`: Chinese version for direct task preparation.
- `doc/README.md` and `doc/README_zh_CN.md`: documentation entry points.
- `doc/source/architecture_review.rst`: points to the task specification rather
  than maintaining a second, partially duplicated prompt.

The specification defines project context and exact base revision, observable
objectives, system-wide ROI prioritization, source-backed pre-change
investigation, implementation constraints, verification evidence levels,
documentation and scientific provenance, traceability, scope control,
acceptance criteria, and delivery records. The architecture review remains the
system-level priority framework; the task template governs individual task
execution; this file remains the current-state handoff.

R1 investigation reconfirmed that `ReadControlInfo(caseDir)` is a production
pipeline with ordered manifest discovery, multi-file parsing, typed database
commit, process coordination, database broadcast, and diagnostic output. A
valid serial fixture must initialize the existing parallel state explicitly;
static-zero fields are not a valid runtime contract. The existing focused
targets compile selected source files and do not yet provide a cohesive target
for the complete reader path. Do not copy the entire executable source list
into a test target. Continue mapping the dependency closure and identify the
smallest production-cohesive integration target.

Current implementation status remains:
- Manifest discovery has component-level coverage, including append and order
  semantics.
- R1 remains Partial: production multi-file reader integration, effective
  `SetupCaseEnvironment()` path application, sequential-case isolation, and
  failure cleanup remain unverified.
- No local build or CTest run was performed in this environment for the latest
  documentation and test changes. The 977/977 result is the user's pre-change
  baseline only.

Next action: inspect exact-HEAD CI for the test-contract change; then map the
production source/dependency closure for `ReadMultiScriptFiles` and
`ReadControlInfo(caseDir)`. Add an integration target only when its source
set and runtime initialization represent a real reusable boundary. Keep the
whole-system priority review active: runtime lifecycle first, mesh topology
and derived-data validity next, then numerical verification and documentation
automation as enabling work.


## 2026-10-30 continuation: ingestion dependency boundary and CI checkpoint

Starting revision: `4b2e5d496ed359c1e3c3645413c4541dc81122f7` on `hexin`. The user reports 977/977 local tests passed before this continuation; this is a pre-change baseline only.

The source-level dependency map and revised test layering are recorded in [the R1 integration investigation](plans/setup-case-environment-integration-test.md#2026-10-30-follow-up-source-level-dependency-graph-checkpoint), committed as `10a21b46b3f071f233dcae854d24dad09c1c91c2`. The key decision is to test the real configuration-ingestion composition first (manifest -> multiple production script parsers -> typed DataBase commit), separately from the broader `ReadControlInfo(caseDir)` orchestration (parallel coordination, broadcast, logging, PIO and database dump). This avoids turning the first test into an improvised miniature of the full executable.

Exact-HEAD CI observed for starting revision `4b2e5d496ed359c1e3c3645413c4541dc81122f7`:
- Sphinx build: **success** — [job](https://github.com/eric2003/OneFLOW/actions/runs/41059549136/job/124075558246).
- Linux OneFLOW Actions: **success** — [job](https://github.com/eric2003/OneFLOW/actions/runs/41059549091/job/124075553594).
- Windows OneFLOW Actions: **in progress** at the time of inspection — [job](https://github.com/eric2003/OneFLOW/actions/runs/41059549210/job/124075554479).

No source code or tests were changed in this continuation; no local build or CTest was run. The documentation-only dependency-map update is now at `10a21b46b3f071f233dcae854d24dad09c1c91c2`. R1 remains **Partial**. Next: implement the real multi-file configuration-ingestion integration target with RAII-managed temporary fixtures, actual DataBase commits, typed assertions, and explicit failure semantics. Then add the broader `ReadControlInfo(caseDir)` orchestration test with a correctly initialized serial/MPI runtime and diagnostic output. Recheck CI on the new exact HEAD before claiming validation.


## 2026-10-30 continuation: production configuration-ingestion composition

Starting `hexin` revision: `4b359903f4ffbf3bf6c595f066d920bf6b26e5b3`. User-reported local baseline: 977/977 tests passed before this change; do not present it as post-change validation.

Implemented on `hexin`:
- Moved `ReadMultiScriptFiles` from `ControlFileReader.cpp` into `ControlScriptParser.cpp`, separating ordered parsing/commit from runtime parallel orchestration and database dumping.
- Changed its input to `const std::vector<std::string>&` and documented ordering plus per-file incremental commit semantics in `ParaFile.h`.
- Added `tests/io/control_configuration_ingestion_test.cpp`: a unique temporary case with a real manifest and two scripts is read through the production manifest reader and multi-file parser, then typed values and later-file replacement are asserted against the actual global DataBase.
- Expanded and renamed the combined CMake target to `control_configuration_test`; existing resolver/adapter tests remain included. The target links the production parser, manifest, resolver, adapter, database and low-level support files, but deliberately excludes `ControlFileReader.cpp` and its process/broadcast/dump dependencies.
- Detailed implementation scope and limitations are recorded in `doc/plans/setup-case-environment-integration-test.md` (commit `8ef8c740fbb27e85d1ea5eb9f83950da8a2097b6`).

Validation status: **Implemented, unverified**. No local compiler/CTest run was available here. The 977/977 result is the pre-change baseline. This test does not cover `ReadControlInfo(caseDir)`, parallel initialization, broadcast, database diagnostic output, `SetupCaseEnvironment`, sequential-case isolation, or failure cleanup. R1 remains **Partial**. Next: inspect exact-HEAD Windows/Linux CI, resolve any CMake/link errors, request and record the user's post-change local test result, then add a separate faithful orchestration test only after its runtime dependency closure is explicit.


## 2026-10-30 continuation: CI regression diagnosis and CMake correction

Starting revision: `326a576e9271ed00efaca30a769d70b87f18ebe5`. User-reported local result remains 977/977 from before the changes; it is not post-change evidence.

Exact-HEAD CI for `326a576e...` exposed the same configuration failure on both platforms, before compilation or tests:
- Windows configure failed at `tests/io/CMakeLists.txt:59`: CMake could not find `control_configuration_test.cpp`, then reported no sources for `control_configuration_test`.
- Linux configure failed for the identical reason.
- Sphinx build passed.

The source inventory confirms the existing resolver/adapter test file is still named `tests/io/control_parameter_resolver_test.cpp`; the target had been renamed but its source path had not. Corrected the CMake source entry to the existing filename in commit `d9180818eed6946130cfdd29a96477ff8c932e7d`. No production C++ behavior was changed by this correction.

At the time this handoff was written, Windows and Linux CI for `d9180818...` were both in progress:
- [Windows job](https://github.com/eric2003/OneFLOW/actions/runs/41063294279/job/124088158630)
- [Linux job](https://github.com/eric2003/OneFLOW/actions/runs/41063294280/job/124088158905)

The new integration test currently verifies the successful two-file manifest -> production parser -> typed global DataBase path. It does not yet verify missing-file failure semantics, `ReadControlInfo(caseDir)`, MPI/serial orchestration, diagnostic output, case path application, sequential-case isolation, or cleanup at the real case runtime boundary. R1 remains **Partial**. Next action: wait for exact-HEAD CI, inspect any subsequent compile/link/test failures, and only then extend the integration behavior. Avoid treating a successful CMake configure as proof that the expanded source target compiles or runs.


## 2026-10-30 continuation: const DataObject lookup compile fix

User-reported MSVC diagnostic while compiling `control_configuration_ingestion_test.cpp`: `DataBase.h:152` could not resolve `GetDataValue<T>(&data)` when `data` was `const DataObject&`.

Root cause verified in the declarations:
- `DataObject*` overload declared `int iElement = 0`.
- `const DataObject*` overload declared `int iElement` with no default, although its definition takes the same two parameters.
- The const-DataBase `GetDataValue` implementation passes only the `const DataObject*`, so overload resolution fails.

Corrected `codes/database/include/DataObject.h` by giving the const overload the same default element index (`= 0`). This preserves const correctness and aligns the two scalar lookup APIs; it does not cast away constness or change runtime behavior. Commit: `e8fa80ec66ed78790e7e9b13b62e22503823cf94`.

The ingestion integration test's `const DataBase&` lookup now exercises this const overload directly. This is useful API-contract coverage: const database access must support scalar element-zero reads just as mutable database access does.

Validation: source-level fix only at this point. No local MSVC build or CTest has been run by the assistant. The 977/977 result remains the pre-change baseline. Next: verify exact-HEAD Windows/Linux CI; inspect any compile/link/runtime failures before expanding R1. R1 remains **Partial**.


## 2026-10-30 continuation: database contracts and documentation

User-reported local validation before this continuation: 978/978 passed. The previous exact HEAD `cc73e947ceb0c9b9b88477f81aa0578012ba827b` had successful Windows, Linux, and Sphinx checks.

Changes made:
- Added API-level documentation to `DataObject.h` for the type-erased value boundary, borrowed raw storage access, and the preconditions of mutable/const `GetDataValue<T>` (non-null object, exact stored type, valid index). The two overloads have matching default index zero.
- Documented `DataEntry`/`DataPara` ownership, borrowed lookup lifetimes, and same-key update behavior in `DataPara.h`.
- Documented `DataBase` ownership of parameter and field registries, and mutable/const named lookup contracts in `DataBase.h`.
- Extended the real configuration-ingestion test to exercise both mutable and const `DataBase` scalar reads.
- Added `doc/source/database_architecture.rst` and linked it from the Sphinx index. It explains ownership, type erasure, parameter replacement, constness, field pointer registration, serialization, and the boundary between ingestion integration and full runtime integration.

Commits in this documentation/API slice: `e1e3a231226635709156b939f686f47661268ffb`, `8504f6a5c2e85efe146823216f6b19d4ed7fef03`, `afe134cc5a5451a94341fb7b05f458aebbcb3999`, `c4b816cbd865223360f685a23ccd96dd4344d993`, `a6abebde9dc3120240ef2edf713b14a2c8ea8bc5`, `eb8590aef26e2ac0a47c9cb4ed85aa5ef5ca48ed`.

Validation status: these changes have not yet been validated by a new local run or exact-HEAD CI. The 978/978 result is the pre-change baseline for this slice. R1 remains Partial. Next inspect exact-HEAD CI, then continue with explicit missing-file/parse-failure semantics and the separate runtime orchestration layer. Avoid implying the database guide proves unsupported behavior; it records current source-backed contracts and clearly separates remaining gaps.


## 2026-10-30 continuation: configuration-ingestion failure semantics

The user reports a pre-change local baseline of 978/978 tests passed. Do not
treat this as validation of the following commits.

Starting `hexin` HEAD was
`3fadcbacafee01081b65bc889218045566dc6a46`. Sphinx passed on that exact
revision; Windows and Linux build/test jobs were still in progress when
inspected.

Source review of `ConfigLoader`, `ControlScriptParser`,
`ConfigDatabaseAdapter`, and `Prj::OpenFile` established the following
current behavior:

- The multi-file reader is sequential and incremental across files.
- Each file is fully parsed into a fresh `ConfigDocument` before its adapter
  commit begins.
- Missing-file failure occurs during file open and leaves prior file commits
  intact.
- A parse-time failure prevents the current file's document from reaching the
  adapter, so declarations parsed earlier in that same file are not committed.
- The adapter's commit is itself not atomic across entries; a later entry
  failure can leave earlier entries from the same document committed.

Added tests for the first two failure boundaries in
`tests/io/control_configuration_ingestion_test.cpp`, committed as
`7c3d488665a35cd5aa293822c64ce9f8b5b91cd3`. The investigation plan was
updated in `doc/plans/setup-case-environment-integration-test.md`, committed
as `df4391a13de37950b24e3889df36eeb7dec5b65b`.

The tests intentionally specify current behavior; they do not introduce
whole-manifest rollback. No local build or CTest was run in this environment.
The user's 978/978 result is the pre-change baseline only. Recheck exact-HEAD
Windows, Linux, and Sphinx CI after the latest documentation commit.

R1 remains Partial. The next larger boundary is still the real
`ReadControlInfo(caseDir)` orchestration, including initialized serial
parallel state, broadcast path, and `log/database.log`; after that, verify
`SetupCaseEnvironment` path application, sequential-case isolation, and
failure cleanup. Do not expand the test target by copying the full executable
source list.


## 2026-10-30 latest continuation: script reader boundary

User-reported local baseline before this slice: 981/981 passed. The starting
`hexin` HEAD was `0e0054d658db092e1cfbc54195ac87d6ef3caf84`; do not treat
that baseline as validation of this slice.

The case-level ingestion test now calls the production
`ReadPrjScript(caseDir)` composition rather than manually invoking manifest
discovery and multi-file parsing. To keep this production entry point
independently linkable, `ReadPrjScript` overloads were moved from
`ControlFileReader.cpp` into the new
`codes/io/src/ControlScriptReader.cpp`. `ControlFileReader.cpp` now retains
the `ReadControlInfo` orchestration boundary (parallel greeting, broadcast,
and database dump). The new source is included by recursive production source
collection and explicitly by `control_configuration_test`.

This separation is an intentional dependency boundary: script ingestion is
distinct from process orchestration. It does not claim full
`ReadControlInfo(caseDir)` coverage. Exact-HEAD Windows/Linux CI and local
CTest are required; no local build was run here. The 981/981 report predates
these edits.

Next: inspect the exact link dependency closure of `ReadControlInfo` and
consider a cohesive production target only if it is useful to both executable
and tests. Preserve real serial/MPI initialization, DataBook broadcast,
logging, PIO, and database dump semantics. R1 remains Partial until
`ReadControlInfo`, effective Case path application, sequential-case
isolation, and failure cleanup are verified.


## 2026-10-30 continuation: production reader orchestration slice

The user reports 981/981 tests passed before this slice. This is a pre-change
baseline only. Starting `hexin` revision:
`144cde3f472fc73d2a8633f84b5eaea874308de8`.

Implemented source changes:
- Corrected the case database dump to use `PIO::OpenCaseFile` rather than
  calling `Prj::OpenCaseFile` directly. The latter performs path resolution
  and file creation but no parallel-rank filtering; the old case dump could
  therefore open/truncate the same output on every MPI rank while its close
  operation was owner-filtered.
- Split explicit path-based PIO operations into
  `codes/io/src/PIOFileAccess.cpp`. The remaining `PIO.cpp` operations are
  task-state-bound convenience entry points. This is a real dependency
  boundary, not a duplicate test implementation.
- Extended `control_configuration_test` to compile the production
  `ControlFileReader.cpp`, `ControlDatabaseIO.cpp`, `PIOFileAccess.cpp`,
  and `LogFile.cpp`.
- Added a serial integration test that initializes the production serial
  parallel state, invokes `ReadControlInfo(caseDir)`, and checks manifest
  ordering/override, typed global DataBase values, `log/database.log`
  contents, and the case-bound diagnostic log.

The new test covers the real serial reader orchestration and the serial
broadcast short circuit. It does not establish MPI collective correctness or
multi-rank file ownership at runtime. Full `SetupCaseEnvironment()`,
effective path application through that entry point, sequential-case
isolation, and exception cleanup remain unverified. R1 remains Partial.

The source-backed dependency and ownership analysis is recorded in
`doc/plans/setup-case-environment-integration-test.md`.

Validation: no local compiler or CTest run was performed here. The user's
981/981 result predates this slice. Inspect exact-HEAD Windows/Linux CI and
request a fresh local test result before claiming the new integration test
passes.


## 2026-10-30 continuation: geometry-core PIO link dependency

The user reported MSVC linker errors from `GridSummaryAdapterTest.exe`:
`OneFLOW_GeometryCore.lib(Zone.obj)` and `GridGroup.obj` referenced
explicit-path PIO methods that were no longer defined in `PIO.cpp`.

Root cause: `tests/geometry/CMakeLists.txt` builds the shared
`OneFLOW_GeometryCore` static library from an explicit source list. After
explicit-path PIO methods were moved to `codes/io/src/PIOFileAccess.cpp`,
this target still listed only `PIO.cpp`. The production executable's
recursive source collection was not affected, and the IO ingestion test
already listed the new source; the geometry test core was the missing
consumer.

Added `PIOFileAccess.cpp` beside `PIO.cpp` in the geometry-core source list.
This fixes the shared library dependency for all geometry tests linking
`OneFLOW_GeometryCore`, rather than patching only `GridSummaryAdapterTest`.

Commit: `eab78cbcc4a0b9ee1fb59a439844a67c6bf7a687`.
No local build or CTest was run here. Reconfigure/rebuild the Windows Debug
geometry test target and run the full suite; this commit is not yet validated
by the user's local test result. Also inspect exact-HEAD CI before treating
the branch as green.


## 2026-10-30 continuation: Case lifecycle integration boundary

The user reports 982/982 local tests passed before this continuation. A fresh
source review mapped the real production path across
`Simulation::RunImpl`, `SimuImp::RunCase`, `SimuContext::SetupCaseEnvironment`,
`TeardownCase`, and process finalization.

Important ownership result: one process context owns parallel/accelerator
initialization and finalization; each case gets a distinct `SimuImp` /
`SimuContext`; case setup and teardown run inside the case cleanup boundary.
A failure in `SetupCaseEnvironment` still reaches `PostProcess` /
`TeardownCase`, after which the outer process lifecycle finalizes the shared
runtime.

The IO path reset contract already has a component test:
`PrjIOPaths.DomainRootsResetWhenBindingNextCase`. It configures project and
domain-specific roots on one case, rebinds a second case, and verifies all
mesh/flow roots resolve under the second case. This proves the `Prj` path
component reset, not the complete `SetupCaseEnvironment` composition.

The remaining R1 gap is therefore not missing isolated tests for each
operation: the real `ReadControlInfo`, IO settings application, and
`Prj` root reset are individually exercised. What remains untested is their
production composition inside `SimuContext::SetupCaseEnvironment`, plus
cleanup after failure at that boundary. The investigation and next extraction
proposal are recorded in
`doc/plans/setup-case-environment-integration-test.md`.

No source code or post-982 test result is claimed in this slice. The latest
hexin revision for this note is recorded by its commit history; check exact
HEAD CI before declaring the preceding PIO source-list fix green.


## 2026-11-01 continuation: shared production runtime and SimuImp integration

The user reports 1096/1096 local tests passed before this slice. This is a
pre-change baseline only.

Starting `hexin` HEAD: `0da2d68deeb7216c174e82630d75dab792103722`.

Implemented on `hexin`:
- `codes/CMakeLists.txt` now builds the recursively collected production
  source set as static `OneFLOW_Runtime`; `OneFlow.cpp` remains the sole
  executable entry point. The runtime target publishes the prior executable
  compile and link requirements so production consumers share one source graph.
- `SimuContext::IsProcessEnvironmentActive()` exposes a read-only process
  lifecycle observation.
- `simu_imp_case_lifecycle_integration_test` links the production runtime and
  invokes real `SimuImp::RunCase()` with a temporary case whose second control
  script is missing. It checks the full Context teardown through seeded task
  and solver-name state, configuration rollback, and process-runtime lifetime.
- Runtime architecture, lifecycle investigation, roadmap, and this handoff
  document the target boundary and acceptance gaps.

Important test contract: after the failed case, task resolution and expanded
solver names must be reset, partial configuration and project-root binding
must be absent, and process runtime must remain active until explicit
`FinalizeEnvironment()`.

No local build or CTest was run here. The 1096/1096 report is pre-change
baseline only. Exact-HEAD Windows/Linux CI and a fresh local build/test are
required because the executable was moved behind a new static-library target.

R1 is improved with real production `RunCase() -> PostProcess() ->
TeardownCase()` integration coverage, but remains Partial. Remaining:
individual teardown-step failure injection in the concrete chain, successful
repeated full SimuImp cases, multi-rank MPI, and CFD numerical regression.
Next action: inspect exact-HEAD CI and fix any build/link regressions before
expanding the integration suite.


Latest validation poll for `db7439f55c9b230004e09640f12f31fa860bcfad`:
- Windows CI queued: https://github.com/eric2003/OneFLOW/actions/runs/41281608533
- Linux CI queued: https://github.com/eric2003/OneFLOW/actions/runs/41281608515
- Earlier CI runs for the preceding C++ commits were still in progress; none had completed at this poll.
- No local CMake configure, build, or CTest was run in this environment. The user's 1096/1096 report remains the pre-change baseline.


## 2026-11-01 continuation: Verification Contract Workbench proposal

The user reports **1125/1125 local tests passed** on the preceding revision. The latest inspected hexin baseline before this documentation slice was the merge commit `74baec02d1e170340224c6e320b02412f184947d`.

The current branch already contains two useful foundations:
- Numerical Contract Workbench: `tools/verification/curve_catalog.py` and `tools/oneflow_gui/verification_curve_dialog.py`, with formula/principle/reference/source/API/documentation/code-sketch/curve data and explicit proves/does-not-prove statements.
- Geometry Contract Workbench: `tools/geometry_workbench` and `tools/oneflow_gui/geometry_contract_dialog.py`, with editable fixture-backed scenarios. Its current stub does not execute the production geometry algorithm; keep that limitation visible.

Added a discussion proposal:
- `doc/plans/verification-contract-workbench.md`
- Indexed from `doc/source/index.rst`.

The proposal recommends a shared contract envelope with domain-specific typed payloads and evidence adapters. Candidate domains include numerical formulas/methods, geometry/mesh, runtime lifecycle, workflows, case regressions, I/O/configuration, and performance. Common navigation is Overview, Principle, References, Architecture, Implementation, Code, Documentation, and Verification, with optional type-specific visualizers such as Formula/Curve, Scenario, Case, Workflow Graph, Mesh View, Benchmark, and Run Artifacts.

Important design decisions proposed (not yet implemented):
- A Case is normally an evidence fixture linked to multiple contracts, not itself the contract.
- Definitions remain versioned and reviewable in the repository; generated run evidence/artifacts live under user data paths.
- Begin read-only and preserve existing source catalogs; do not immediately migrate everything to JSON/YAML.
- Keep catalog/schema/evidence resolution Qt-independent.
- Do not add generic dynamic plugin execution. Use typed, allow-listed adapters that invoke existing test and case systems.
- Preserve the explicit evidence scope: a passing unit test does not establish full physical validation.
- Keep geometry fixture-backed results distinct from production algorithm results.

Industry practices to borrow selectively: requirements traceability, V&V separation, reproducible run manifests, layered testing with production-path integration, ADRs, schema evolution, and explicit claim limits.

Status: **Proposal committed; no runtime or Workbench implementation changes in this slice.** User-reported 1125/1125 is the baseline, not validation of new executable code. Next step is to discuss the five open design questions in the proposal, then inventory current numerical and geometry metadata before Stage 1. Do not start a broad rewrite.


## 2026-11-01 decision: incremental Contract Workbench evolution

The user agrees with a gradual migration: **keep the existing workbenches running first, then abstract and refactor incrementally**. Do not pause working functionality for a broad rewrite.

Working rules:
- Preserve the current Numerical Contract Workbench and Geometry Contract Workbench entry points and behavior while adding shared capabilities.
- Start with read-only compatibility/adaptation; do not migrate their source data formats prematurely.
- Extract shared abstractions only after at least two real consumers demonstrate the same stable need.
- Keep execution in authoritative existing systems (GoogleTest/CTest, geometry native tests, CFD regressions, workflow runners). The Workbench organizes and presents evidence; it must not duplicate solver algorithms or become a competing test engine.
- Keep the shared catalog/validation core Qt-independent and test it headlessly.
- For each migration step, maintain old-path/new-path equivalence tests and a runnable fallback until the new path is established.
- No large-bang GUI, catalog, or runner rewrite.

User-reported validation baseline for this continuation: **1125/1125 local tests passed**. This validates the preceding user-tested revision, not any future implementation work.

Recommended first implementation slice after the design review:
1. Inventory current Numerical entries and Geometry scenario metadata.
2. Define the smallest read-only common view/protocol needed to list a contract's stable ID, title, domain/kind, overview/principle, references, implementation/documentation links, and evidence links.
3. Adapt one existing numerical entry and one geometry scenario without changing either source format or current UI.
4. Add headless tests for adaptation, missing references, and explicit fixture-vs-production evidence labeling.
5. Only then add a shared read-only browser surface or common navigation shell.


## 2026-11-01 implementation: first read-only Contract adapters

The user confirmed the migration strategy: keep existing workbenches working,
then abstract and refactor incrementally. The first executable-code slice now
exists on `hexin`:

- `tools/verification/contract_catalog.py`
  - Immutable, Qt-independent `ContractSummary` and `ContractEvidence`.
  - `adapt_theory_entry` adapts the current numerical dataclass without
    changing `curve_catalog.py`.
  - `adapt_geometry_scenario` adapts an already loaded geometry scenario
    without changing its JSON format or running either the fixture pipeline
    or production C++.
  - `validate_contract_summaries` checks duplicate IDs, path containment,
    and existence of source, documentation, theory CSV, and linked-test
    source files.
- `tools/verification/test_contract_catalog.py`
  - Covers both adapters, evidence-status semantics, invalid/duplicate IDs,
    missing and escaping paths, and repository references.
- `tools/verification/README.md` and
  `doc/plans/verification-contract-workbench.md` describe the compatibility
  layer and focused validation commands.

Evidence semantics are deliberately conservative:
- numerical API names and theory CSV are `declared_not_executed`;
- a geometry fixture is `fixture_available`, not production verification;
- linked GoogleTest declarations are `linked_not_executed`;
- a non-null scenario result is `result_recorded_unverified` until a
  trustworthy runner/provenance model exists.

No existing GUI entry point or source data format has changed. This slice
does not add a generic runner, JSON/YAML migration, or shared GUI shell.

Assistant has not run the new Python tests or full CTest in the user's
Windows/MSVC environment. The user's 1125/1125 result remains the
pre-implementation baseline. First run:

```powershell
python -m unittest tools.verification.test_contract_catalog
python -m unittest tools.verification.test_curve_catalog
```

Then run the full existing suite. Resolve any failures before adding more
consumers. Next likely increment: make the adapter inventory all current
numerical entries and geometry scenarios, report invalid/stale links, and
only then expose a read-only combined browser view. Do not change the two
existing GUI flows yet.

## 2026-11-01 implementation: combined contract inventory

The next incremental slice extends the read-only compatibility model:

- `load_contract_catalog(repository_root)` enumerates all current Numerical
  `TheoryEntry` records and all Geometry scenario JSON files through their
  existing source loaders.
- Canonical IDs are namespaced as `numerical/<source-id>` and
  `geometry/<source-id>`; `source_id` retains the original identifier.
- Geometry fixture evidence records the repository-relative JSON path.
- Validation checks duplicate IDs, repository-relative path containment,
  source/document/data file existence, and whether each linked GoogleTest
  declaration is actually present in its declared C++ source.
- `tools/oneflow_gui/test_test_inventory.py` explicitly verifies that the
  new `tools/verification/test_contract_catalog.py` methods are discovered
  by the existing Workbench test inventory.

This remains read-only discovery. It does not invoke GoogleTest, CTest,
geometry production algorithms, or persist execution results. The test
inventory's discovery of a Python test is not evidence that the test passed.

Focused commands:

```powershell
python -m unittest tools.verification.test_contract_catalog
python -m unittest tools.oneflow_gui.test_test_inventory
python -m unittest tools.verification.test_curve_catalog
```

The assistant has not run these commands in the user's Windows/MSVC
environment. The last user-reported 1125/1125 remains a pre-change baseline.
Do not add GUI integration until these checks pass. The next design step after
validation is to inspect the generated catalog and decide which traceability
fields are genuinely missing before introducing Runtime/Workflow contracts.

The combined loader also reuses the existing Numerical catalog's pure-API
reference validator. It checks that declared API references resolve to
repository symbols, in addition to validating source/data/document paths and
linked Geometry GoogleTest declarations. This reuses established validation
logic rather than creating a second parser for the Numerical catalog.


## 2026-11-01 implementation: Contract Workbench tests in GUI Unit Tests

The user reported 1125/1125 local CTest tests passed before this slice and
confirmed these focused Python modules pass locally:
- `tools.verification.test_curve_catalog`: 5 tests.
- `tools.verification.test_contract_catalog`: 9 tests.
- `tools.oneflow_gui.test_test_inventory`: 6 tests.

Reviewed the current `hexin` GUI test architecture:
- `test_inventory.py` performs read-only source discovery of GoogleTest
  declarations under `tests/` and Python unittest methods under `tools/`.
  It can attach Geometry scenario input/expected metadata only when the
  scenario explicitly links a test. Its result field remains unset.
- `repository_test_inventory_dialog.py` presents that catalog; it does not
  execute tests.
- `test_runner.py` owns an explicit Python command registry plus the CTest
  command. `test_dialog.py` executes these commands sequentially.
- The GUI's Unit Tests action therefore differs from the Repository Test
  Inventory: one executes registered commands; the other reports declarations.

Changes on `hexin`:
- Added the three passing modules to `_PYTHON_UNIT_TESTS` in
  `tools/oneflow_gui/test_runner.py`.
- Extended `tools/oneflow_gui/test_test_runner.py` to assert that all three
  module paths remain registered in the GUI unit-test command list.
- Added `doc/plans/testing-architecture.md` and linked it from the Workbench
  README. The document describes test-layer ownership, each module's purpose,
  discovery-versus-execution semantics, evidence limitations, and industry
  practices.

Design decision:
- Unify the GUI entry point and command/result presentation, not the underlying
  test frameworks. CTest remains authoritative for CMake-registered native
  tests; Python unittest remains authoritative for Python tests; the CFD
  regression harness and Development Test case runner retain their own
  execution semantics.
- Keep the source inventory broader than the execution allow-list. Do not
  automatically execute every AST-discovered `test_*.py`: commands may need
  specific arguments, dependencies, or working directories.
- Contract fixture metadata and linked-test declarations remain unexecuted
  until an actual runner reports a result.

Commits:
- `22201a1a23cf7eac91689fe87bd3e2f549b3b894`: register the three modules.
- `c3d5c50709770bde07ea6aefcf2bd4fd07a0a0f0`: test the command registry.
- `1ff595bf6d338231511bccd68a133bfc966b930d`: add testing architecture doc.

No tests were executed by the assistant in the user's Windows environment.
The user's 1125/1125 result and 5/9/6 focused Python results are pre-change
validation for the code/documentation slice. Next validation:
```powershell
python -m unittest tools.oneflow_gui.test_test_runner
python -m unittest tools.verification.test_curve_catalog
python -m unittest tools.verification.test_contract_catalog
python -m unittest tools.oneflow_gui.test_test_inventory
ctest --test-dir build -C Debug --output-on-failure
```


## 2026-11-01 follow-up: split CTest and Python Unit Tests in GUI

User requested separate GUI selections for CTest and Python tests.

Changes on `hexin`:
- Added `build_ctest_test_commands()` and
  `build_python_test_commands()` to `tools/oneflow_gui/test_runner.py`.
  The legacy combined `build_unit_test_commands()` remains for callers that
  need the combined sequence.
- Replaced the single "Unit Tests (CTest and Python tests)" checkbox in
  `test_dialog.py` with independent "C++ Unit Tests (CTest)" and
  "Python Unit Tests (unittest)" checkboxes. Both default to checked, but
  either can now run independently.
- If CTest is selected but `build/CTestTestfile.cmake` is absent, the GUI
  reports that the CMake build tree must be configured; it does not silently
  skip the selected test group.
- Added command-builder tests for independent command construction and
  unconfigured CTest handling.
- Updated GUI and testing architecture documentation.

Commits:
- `87c393d41c3b25f4bfd348f8d4640774b3204817`: add independent command builders.
- `d15565aab3e5e47021bb10b33cd194f38a50bfeb`: split GUI test selections.
- `7d178caa3f7c7906b11471d1e3a4ed41a3ad8348`: remove combined builder from GUI path.
- `f6f280e7e1a74e0c9f3f6e28299a3d12e282aa44`: test independent command construction.

Assistant has not run tests in the user's Windows environment. Next validation:
```powershell
python -m unittest tools.oneflow_gui.test_test_runner
python -m unittest tools.oneflow_gui.test_test_inventory
python -m unittest tools.verification.test_curve_catalog
python -m unittest tools.verification.test_contract_catalog
```
Then open the GUI Tests dialog and run CTest and Python Unit Tests separately.

## Python test inventory and execution-registry audit (2026-11-01)

The Workbench Python test command registry was expanded from 24 to 41
commands. Newly registered coverage includes the mesh-layout editor document
and rectangle compiler, Geometry Workbench analytic/scenario tests, background
task services, project/template inventory, SA closure curves, and the standalone
Tkinter runner test.

tools/oneflow_gui/test_inventory.py now discovers Python unittest methods
under tools/, tests/, and gui/python/tkinter/. Its
audit_python_test_registry() compares discovered test source files with
explicit module, discovery, and script commands in
tools/oneflow_gui/test_runner.py. The audit reports discovered, registered,
excluded, unregistered, and stale-exclusion paths. Every discovered module
must be registered or have a documented exclusion.

Current intentional exclusions:
- Qt widget/queue and dialog tests pending a dedicated headless GUI test policy.
- LaTeX preview tests pending the optional Matplotlib rendering environment.

The audit is exercised by
tools.oneflow_gui.test_test_inventory, including a negative test proving
that an unregistered unittest module is detected. Example-local/tutorial
verification scripts are outside the generic unit-test roots; they retain
their own case-specific execution contracts.

Validation still required on the user's Windows checkout:

    python -m unittest tools.oneflow_gui.test_test_runner
    python -m unittest tools.oneflow_gui.test_test_inventory
    python -m unittest tools.oneflow_gui.test_test_runner tools.oneflow_gui.test_test_inventory

Run the full Python Unit Tests group in the GUI after pulling the changes.
The assistant did not execute the updated suite in the user's Windows
environment.

Follow-up cleanup from the user's GUI test output:
- workflow_process.run_process() now closes the captured subprocess stdout
  stream when the reader thread reaches EOF, preventing an unclosed-file
  ResourceWarning. The Windows process test now uses io.StringIO so the
  mocked stream follows the same closeable-file contract.

## Background task stop-state correction (2026-11-01)

The user's Windows run exposed a failure in
test_background_task_host.test_host_stops_active_process_when_stop_file_appears:
the task ended as failed instead of stopped. The host now latches the stop
request as soon as the stop file is observed, performs process-tree termination
and fallback cleanup, and finalizes as stopped when the requested process has
been reaped. The test explicitly verifies that its stop-request thread delivered
the request before asserting the persisted terminal state.

Validation required on Windows:

    python -m unittest tools.oneflow_gui.test_background_task_host
    python -m unittest tools.oneflow_gui.test_background_task_manager


Follow-up hardening on the same slice:
- Contract JSON sections are shape-checked before adaptation.
- Declared IDs must use their domain as the namespace.
- Declared evidence statuses are validated against explicit semantics
  (`linked_not_executed`, `fixture_available`, `passed`, `failed`,
  `blocked`, `stale`, etc.); arbitrary labels such as `green` are rejected.
- Test fixtures were updated for the common `architecture` and `status`
  fields.
- Latest focused-test source commits:
  `69e362b50a28efa1c95a783f0ba300c6006082b8`,
  `b19040a57bb6338a6f95ef12bf7379ddb48353ce`.
- Latest catalog validation commits:
  `e0dc99728de0cad69244005c01a80b11de7a3033`,
  `2bb1d6af8ec65156245c346b211726c8cf661bf8`,
  `70230b22886f4ce9fbf69af4d29bea9a9f33e48c`,
  `7f23f44ab972506bb1e3c8da6b64a8611edb19de`.

The user's 1126/1126 result was reported before these Contract Workbench
changes. Re-run the focused catalog tests and full suite; do not mark this
slice verified until then.


Final evidence-link correction: declared `google_test` records are now
validated by the same source-declaration checker as Geometry
`linked_test` records. This ensures the Runtime contract's three test links
are not merely path-checked; each named GoogleTest suite/case must be present
in its declared C++ source. Commit:
`1f28288bfa88cf87d1d86bbb9eb76d346bcfab8a`.


Case-contract follow-up hardening after the initial handoff entry:
- `CaseContractData` now explicitly declares residual output channels.
- The catalog requires residual outputs to be a subset of the compared outputs.
- It validates that the output manifest and contract output lists match exactly.
- It checks that each residual output exists in the canonical residual database,
  has `iter` and `sub-iter` columns, and has exactly `max_steps` rows.
- It rejects unsupported Case kinds and malformed residual database roots.
- Tests include tolerance drift and residual-history length drift.
- Documentation records these additional invariants.

Latest focused-slice commits include:
`3c8882f4e3afd2f685f2ef06dab539bd98504098`,
`fe471a278ec9d08e9812c5302719ff9b8053268a`,
`be6ea51bb09e0fc4cd91e1401e07d0394de88ae1`,
`8674ede4725e12e0cd769244fb455ab04295e9fa`,
`f5657d890e718cc3ecc63f7a625f7fed3e239056`,
`bd8bf7ef0cb1fb881307c10f7f00630906a05c3b`.

The exact latest hexin HEAD at this handoff is
`bd8bf7ef0cb1fb881307c10f7f00630906a05c3b`. No CI status was reported for
this direct branch commit and no local Windows test run is claimed.


## 2026-11-01: visible cross-domain Contract Workbench entry

The user explicitly requested that the Contract Workbench be visible in the
desktop interface, with distinct entry points for Formula/Knowledge, Geometry,
Case, and Runtime. This request supersedes the earlier temporary instruction
to avoid expanding `tools/oneflow_gui` for this specific feature. The repository
root `gui/` remains out of scope.

Implemented on `hexin`:
- Added `tools/oneflow_gui/contract_workbench_dialog.py`.
- Added a domain hub with Formula & Knowledge, Geometry, Case, and Runtime
  buttons.
- Formula & Knowledge routes to the existing
  `VerificationCurveDialog` (the existing title already identifies it as
  Contract Workbench — Formula knowledge).
- Geometry routes to the existing `GeometryContractDialog`.
- Case and Runtime open a shared read-only `ContractDomainDialog` driven by
  `load_contract_catalog()`.
- Case view adds a typed Case Setup tab for suite, solver, physics, mesh/input,
  compared outputs, residual policy, and analytical reference paths.
- Runtime uses the shared Overview, Principle, Architecture,
  Implementation/Code, Documentation/References, and Verification tabs.
- Replaced the two separate home-page Geometry Contracts and Verification
  Curves buttons with one **Contract Workbench...** entry.
- Added `test_contract_workbench_dialog.py` for domain buttons and Case/Runtime
  navigation, and registered it as the 42nd Python command in
  `test_runner.py`.
- Updated `tools/oneflow_gui/README.md` and the Chinese
  `doc/source/workbench_user_guide_zh.rst`.

The new UI is a presentation layer only. It does not execute C++ tests, run CFD,
or upgrade evidence statuses. Existing Formula and Geometry domain tools remain
in place; this is an incremental integration rather than a rewrite. Workflow
and Performance buttons are intentionally deferred until their typed contracts
have grounded repository content.

Focused validation to run locally:
```powershell
python -m unittest tools.oneflow_gui.test_contract_workbench_dialog
python -m unittest tools.verification.test_contract_catalog
python -m unittest tools.verification.test_curve_catalog
```
Then run all registered Python commands and the full CTest suite. The previous
1126/1126 and 41-command pass results predate this UI integration; no local
Windows test result is claimed for these commits.

The new UI test is registered as the 42nd Python test command. The latest
`hexin` HEAD for this handoff is
`246b1f72ff91af7c7352ec73e5b63309ba4ee876`.


## 2026-11-01: Case Contract artifact navigation

Extended the visible Case Contract view in
`tools/oneflow_gui/contract_workbench_dialog.py`:
- **Open case directory...** opens the declared repository case folder.
- **Open suite and references...** opens the declared suite manifest, residual
  baseline, and retained Blasius reference files.
- Both actions resolve declared paths under the repository root. Directory
  navigation rejects a path that escapes the repository; file navigation uses
  the existing containment-checked `_open_paths` helper.
- These are read-only navigation actions; they do not run CFD or mutate
  reference data.

Extended `test_contract_workbench_dialog.py` to mock the desktop URL-opening
boundary and assert the exact case directory and all five suite/reference
targets. The test remains in the explicit Python command registry (42 commands).
Updated the English Workbench README, Chinese user guide, and Contract Workbench
design plan.

Focused validation:
```powershell
python -m unittest tools.oneflow_gui.test_contract_workbench_dialog
python -m unittest tools.verification.test_contract_catalog
```
The latest user-reported 1126/1126 baseline predates these changes. No local
Windows test execution is claimed.


The desktop-open tests replace the module-level `QDesktopServices` dependency
with a mock object rather than patching a method on the PySide C++ extension
type. This keeps the test independent of Qt's type mutability and the host
desktop shell.

Latest `hexin` HEAD after this slice:
`925a8c779b6944eecbf4f880af46730bc5e2eafb`.


## 2026-11-02 continuation: repository CFD case setup-only integration

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

CMake defines ONEFLOW_SOURCE_DIR for the production runtime integration
test. The contract declaration links the exact GoogleTest name as
linked_not_executed. Static catalog resolution, CMake target declaration,
CTest discovery, actual execution, and pass outcome remain distinct evidence
states. The assistant has not run the user's Windows build or CTest.


## 2026-11-03 continuation: formula preview layout and CFL ramp contract

Based on the current `hexin` source after the 2026-11-02 integration slice:

- Moved the Parameter Formula tab's preview below the symbol definitions, LaTeX source, and formula actions. The GUI test asserts the widget order so this remains an intentional layout contract.
- Extracted `CflAtIteration` as a pure policy from `TimeStep::CalcCfl`. The caller still selects `outerSteps`, or `innerSteps` when `dualtime == 1`.
- Added six focused GoogleTests for ramp start, interpolation, boundary, post-ramp behavior, descending schedules, and a zero-length ramp at iteration zero.
- Added the `parameter/cfl-ramp` Contract Workbench entry and numerical documentation. The source trace distinguishes CFL scheduling from `global_dt`: the local branch uses CFL-scaled cell spectral radii; the global branch uses `ctrl.pdt`; the local-global branch reduces local steps to a common minimum. The supported HIP fast path is restricted to local inviscid stepping.
- Added the CFL ramp to the Run Control group and updated the catalog and GUI search tests.

Validation status: test declarations and source links are checked into the catalog, but they remain `linked_not_executed`. No local Windows/MSVC build, CTest, or Python GUI/catalog test run is claimed for this slice. The latest queried GitHub combined status and PR-associated workflow list did not yet expose checks for the exact final commit. Re-run the focused policy test, Contract Workbench/catalog tests, and full Debug CTest on the exact `hexin` revision before promoting execution evidence.
