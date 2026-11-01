# Accelerator Runtime Shutdown Order

Status: Initial ordering fix validated by user at 1124/1124; retryable-finalization change implemented, pending validation
Baseline for current slice: hexin HEAD ba067226fcfefad896648d24bd7daa75b49b3a20
Scope:
- codes/main/src/SimuContextEnv.cpp
- tests/main/simu_imp_case_lifecycle_integration_test.cpp

## 1. System-level question

OneFLOW has three related but distinct lifetimes:

1. Process runtime: MPI/parallel environment and AccelRuntime.
2. Case context: one SimuContext and its per-case backend-state registry.
3. Backend allocations: device buffers owned by backend state objects, plus shared backend objects such as HipFluxBackend.

The required ordering is:

    backend allocations and shared backend objects
        are destroyed before
    AccelRuntime::Finalize()
        which is completed before
    process-level HXFinalize()

For a normal case, SimuContext::TeardownCase() already clears the registry and
releases the shared HIP flux backend before process finalization. However,
FinalizeEnvironment() is also a public process-level boundary and may be called
directly. Before this change it finalized AccelRuntime without independently
clearing its context-owned registry. That allowed a caller that skipped
TeardownCase() to leave backend state objects alive after the allocator runtime
had been finalized.

This is an ownership-ordering gap, not a numerical solver defect.

## 2. Contract

- SimuContext owns its EulerDomainStateRegistry.
- The registry uniquely owns each EulerDomainState through std::unique_ptr.
- Backend state objects may own device allocations; those allocations require
  the matching backend runtime to remain initialized while they are released.
- TeardownCase() remains the primary per-case cleanup boundary and continues
  to release the registry before solver teardown and process finalization.
- FinalizeEnvironment() is also defensive: before finalizing AccelRuntime, it
  clears any remaining context-owned backend states and releases the shared HIP
  flux backend when HIP is enabled.
- Context-state and shared-backend releases are independent prerequisites and
  are all attempted through ExecuteAllCleanupSteps(). Accelerator runtime
  finalization is a dependent boundary: if it fails, HXFinalize() must not run.
  The lifecycle phase is marked Finalized only after both runtime layers finish.
- This is not a replacement for TeardownCase(). Process finalization does not
  promise to reset every solver/global case object; normal case execution must
  still run the full case teardown path.

The intended partial order is:

    SimuContext::TeardownCase (normal case path)
        -> ClearAccelStates
        -> HipFluxBackend::ReleaseShared
        -> other case-owned solver resets

    SimuContext::FinalizeEnvironment
        -> ClearAccelStates (defensive; idempotent)
        -> HipFluxBackend::ReleaseShared (HIP; idempotent)
        -> FinalizeAccelRuntime (must succeed)
        -> HXFinalize

    If either pre-runtime release fails, all pre-runtime releases are attempted,
    but accelerator finalization is not entered. If accelerator finalization
    fails, HXFinalize is not entered. Both cases leave the process lifecycle
    marker Active so the owner can diagnose and retry.

## 3. Implementation

SimuContextEnv.cpp now separates independent cleanup from dependent runtime
shutdown. It attempts ClearAccelStates() and, in HIP builds,
HipFluxBackend::ReleaseShared() as one pre-runtime cleanup group. If that group
reports a failure, FinalizeAccelRuntime() is not called. After successful
pre-runtime cleanup, FinalizeAccelRuntime() runs alone; HXFinalize() follows
only if accelerator finalization succeeds. The lifecycle marker is committed
only after both calls succeed.

This deliberately keeps ownership in SimuContext and the existing shared
backend owner. It does not add another global manager or make AccelRuntime
depend on solver-specific types. The dependency direction remains:

    main / SimuContext -> accelerator runtime and backend interfaces

The accelerator runtime remains unaware of EulerDomainStateRegistry.

## 4. Runtime finalization failure semantics

AccelRuntime::Finalize() calls the backend before resetting its owning
unique_ptr or clearing the initialized flag. Therefore a thrown backend
Finalize() leaves the backend object and runtime state published. This is the
basis for retry: the backend interface now explicitly requires a failed
Finalize() to leave the backend valid and safe to retry. CPU Finalize() is a
no-op; HIP Finalize() synchronizes the device and does not discard backend
state, so a synchronization error does not trigger premature host-side
runtime destruction.

A focused AccelRuntime test injects a backend whose first Finalize() throws
before changing state. It verifies that IsInitialized() remains true, Backend()
remains accessible, and a second Finalize() succeeds after the injected fault
is disabled. This pins the runtime-level retry contract; it does not simulate
an actual HIP driver/device failure.

## 5. Test evidence

The production-linked test
SimuImpCaseLifecycleIntegration.ProcessFinalizationReleasesContextBackendStatesBeforeRuntime
creates a real SimuImp, initializes its process runtime, inserts a test
EulerDomainState into the real context-owned registry, and finalizes the
process without first calling TeardownCase().

The test state destructor records whether AccelRuntime was still initialized.
Assertions require:
- the state destructor observed an initialized runtime;
- the registry is empty after finalization;
- AccelRuntime is no longer initialized;
- the process lifecycle is marked inactive.

This is an integration-level ordering test, not merely a mock of the cleanup
helper. The test fixture must initialize the production project environment with
Prj::Init() and construct SimuImp from an absolute temporary case directory.
A first run exposed that SimuImp resolves its case path in the constructor:
the placeholder relative path failed before process setup because
Prj::current_dir was not initialized. The test now uses the existing
TemporaryCase fixture and the production Prj initialization path.

The HIP-specific shared backend release branch is compile-time guarded and
must also be covered by the HIP build/CI configuration when available.

## 6. Industry practice and rationale

This follows the C++ Core Guidelines:
- R.1: manage resources automatically with resource handles and RAII;
- C.31: release resources acquired by a class;
- C.36/C.37: destructors must not fail and should be noexcept;
- E.4/E.7: design error handling around invariants and state preconditions.

Reference: https://isocpp.github.io/CppCoreGuidelines/CppCoreGuidelines

RAII handles the lifetime of an individual object. It does not by itself prove
that a device allocator or runtime outlives every object that depends on it.
That cross-object dependency must be expressed by owner ordering and verified
at the process/case boundary. This change makes that ordering explicit without
introducing shared ownership or a global shutdown coordinator.

## 7. Validation status and limits

- User-reported validation of the preceding fixture-correction revision:
  **1124/1124 tests passed**. This is the baseline before the current
  retryable-finalization changes.
- Current changes require a fresh VS2022 Debug build, focused
  AccelRuntimeTest, lifecycle integration test, and full Debug CTest on the
  exact updated hexin HEAD.
- Exact-HEAD Windows/Linux CI and documentation CI must be checked after the
  final documentation commit.
- HIP-enabled compilation and hardware execution are separate evidence. The
  CPU integration test proves the context registry ordering but cannot prove
  device-specific deallocation behavior.
- Multi-rank MPI and teardown failures injected into every concrete solver
  release remain separate validation gaps.
