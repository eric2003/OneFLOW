# Accelerator Runtime Shutdown Order

Status: Implemented; local validation reported by the user is pending for this change set
Baseline: hexin HEAD 6e550260578a7ad3490793791f9a26feb1e670ae
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
- Cleanup steps are ordered, independent, and all attempted through
  ExecuteAllCleanupSteps(). The lifecycle phase is marked Finalized only after
  the cleanup sequence completes without a reported failure.
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
        -> FinalizeAccelRuntime
        -> HXFinalize

## 3. Implementation

SimuContextEnv.cpp adds ClearAccelStates() to the process-finalization cleanup
sequence before FinalizeAccelRuntime(). Under ONEFLOW_ENABLE_HIP it also calls
HipFluxBackend::ReleaseShared() before runtime finalization. ReleaseShared()
resets the shared unique_ptr and is safe to call after the normal case teardown
has already released it.

This deliberately keeps ownership in SimuContext and the existing shared
backend owner. It does not add another global manager or make AccelRuntime
depend on solver-specific types. The dependency direction remains:

    main / SimuContext -> accelerator runtime and backend interfaces

The accelerator runtime remains unaware of EulerDomainStateRegistry.

## 4. Test evidence

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
helper. The HIP-specific shared backend release branch is compile-time guarded
and must also be covered by the HIP build/CI configuration when available.

## 5. Industry practice and rationale

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

## 6. Validation status and limits

- User-reported pre-change baseline: 1123/1123 tests passed.
- This change set still requires local VS2022 Debug build, focused integration
  test, and full CTest on the exact updated hexin HEAD.
- CI must be checked for the exact final commit, including Windows and Linux.
- HIP-enabled compilation and hardware execution are separate evidence. The
  CPU integration test proves the context registry ordering but cannot prove
  device-specific deallocation behavior.
- Multi-rank MPI and teardown failures injected into every concrete solver
  release remain separate validation gaps.
