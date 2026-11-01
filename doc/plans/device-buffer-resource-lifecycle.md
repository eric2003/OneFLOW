# DeviceBuffer Resource Lifecycle

Status: Implemented, awaiting local validation
Baseline: hexin HEAD d9e25a7b671e6104c0d116ae5c561e5d12ba09d0
Scope: codes/accel/include/DeviceBuffer.h

## 1. Problem and system-level impact

DeviceBuffer<T> is the ownership boundary for backend allocations used by accelerator
kernels. It is deliberately grow-only: different solver stages can request different
extents, and retaining the largest allocation avoids repeated device allocation/free
operations in the time loop.

The previous growth path allocated a replacement, copied the existing contents, and then
released the old allocation using a raw pointer. If the copy failed, the replacement had no
owner and leaked. The destructor and noexcept move assignment also called a potentially
throwing backend deallocator through Reset(). A HIP error could therefore terminate the
process during stack unwinding or a move assignment.

This is a resource-lifecycle defect, not a numerical algorithm change. The desired result is
that allocation ownership remains explicit even when backend operations fail.

## 2. Contract

- A buffer is non-copyable and move-only. Exactly one buffer object owns its allocation.
- Resize(n) is grow-only; requests not exceeding the current capacity are no-ops.
- Resize rejects byte-count overflow before allocating.
- CopyToHost rejects requests larger than the allocated capacity before calling the backend.
- A size whose byte count would overflow size_t is rejected before allocation.
- If replacement allocation or device-to-device copy fails, the existing pointer and
  capacity remain unchanged. The temporary allocation is released on a best-effort basis.
- After replacement allocation and copy succeed, ownership is committed to the new
  allocation before releasing the previous allocation. If releasing the previous
  allocation reports an error, the exception propagates but the buffer still owns the
  replacement and remains internally consistent.
- Explicit Reset() clears the object's ownership state before asking the backend to
  release the previous allocation. A release error is propagated; the object remains empty.
- Destruction and noexcept move assignment never propagate backend release exceptions.
  Cleanup in those paths is best-effort.
- A buffer must be reset or destroyed while its owning accelerator runtime is initialized.
  If the runtime has already been finalized, the backend is no longer available to release
  the allocation. The non-throwing destructor avoids termination but cannot recover a
  resource after its allocator has gone away.

The last point is an ordering constraint, not a claim that buffers may safely outlive the
runtime. Long-lived owners must be torn down before FinalizeAccelRuntime().

## 3. Implementation and exception guarantees

The replacement sequence is:

1. Validate count * sizeof(T) for overflow.
2. Allocate a candidate buffer and immediately place it under
   DeviceBufferDetail::AllocationGuard.
3. Copy the old contents into the candidate.
4. Publish the candidate pointer and capacity as the buffer's new owned state.
5. Release the old allocation.

The guard provides scope-bound cleanup for the candidate if step 3 fails. Its destructor is
non-throwing because cleanup during exception propagation must not replace the original
failure. The commit-before-release ordering means a reported old-allocation release error
cannot leave the object pointing at an allocation whose ownership has already been
transferred.

Explicit cleanup and destructor cleanup intentionally have different error policies:
callers may use Reset() when they need to observe a backend error; destruction and
noexcept move assignment use a private non-throwing path.

This follows the C++ Core Guidelines resource-management guidance: [R.1, manage resources
with RAII](https://isocpp.github.io/CppCoreGuidelines/CppCoreGuidelines#r1-manage-resources-automatically-using-resource-handles-and-raii-resource-acquisition-is-initialization),
[C.36, destructors must not fail](https://isocpp.github.io/CppCoreGuidelines/CppCoreGuidelines#c36-a-destructor-must-not-fail),
and [C.37, make destructors noexcept](https://isocpp.github.io/CppCoreGuidelines/CppCoreGuidelines#c37-make-destructors-noexcept).
These guidelines establish the ownership/error-handling practice; they do not prove the
behavior of a particular HIP runtime or device.

## 4. Code and test evidence

Implementation:
- codes/accel/include/DeviceBuffer.h
- codes/accel/include/AccelBackend.h defines the backend allocation interface.
- codes/accel/src/AccelRuntime.cpp owns backend selection and runtime lifetime.

Focused tests:
- tests/accel/accel_runtime_test.cpp
  - DeviceBufferTest.FailedGrowthReleasesCandidateAndPreservesCurrentBuffer
    injects a device-copy failure and checks candidate cleanup plus preservation of the
    original pointer and capacity. The same test also verifies byte-count overflow rejection
    and that CopyToHost cannot read beyond the allocated capacity.
  - DeviceBufferTest.DestructorDoesNotPropagateBackendReleaseFailure
    injects a release diagnostic after the fake backend frees the allocation and verifies
    that destruction does not propagate the exception.
  - DeviceBufferTest.ReleaseFailureKeepsCommittedBufferStateConsistent
    verifies that a failed old-allocation release leaves the replacement owned, and that
    explicit Reset() clears the object even when release reports an error.
  - DeviceBufferTest.MoveAssignmentDoesNotPropagateBackendReleaseFailure
    verifies that the noexcept move assignment remains non-throwing when releasing its
    previous allocation reports an error.

The tests use a tracking CPU backend so failure paths are deterministic and do not require
a physical GPU. They verify host-side ownership and exception contracts; they do not prove
HIP allocator behavior, device-memory correctness, concurrency safety, or solver-level
performance.

## 5. Validation and remaining work

- User-reported pre-change baseline: 1119/1119 tests passed.
- Implementation and focused tests have been committed to hexin; local MSVC/CTest
  validation of these new changes is still required.
- Run the focused accel_runtime_test, then the complete Debug CTest suite.
- Run the representative CPU regression suite; no numerical formula or kernel was changed,
  so a DCU performance claim is out of scope for this slice.
- Follow-up design question: make the runtime/buffer destruction order explicit at the
  owning solver-state boundary, especially for static or process-lifetime state. Do not
  silently claim safe cleanup after runtime finalization.
