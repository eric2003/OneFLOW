# OneFLOW project handoff

**Date:** 2026-09-03
**Repository state:** source commit `60681aaa` was clean and synchronized with `origin/master`; this handoff document is pending
**Latest commit:** `60681aaa` — `Add Euler backend GoogleTest coverage`

## 1. Current conclusion

The project work covers one-dimensional Euler performance, CPU/GPU
correctness, unified residual regression, and a maintainable Git workflow.

The previous GPU/MPI implementation is merged upstream. The latest commit adds
GoogleTest coverage for the Euler backend without changing the production
solver or residual path. CPU tests are validated locally, and the optional HIP
target has now been compiled and run on a real Kunshan DCU node with 6/6 tests
passing.

Do not confuse these states:

- existing Euler CPU/HIP functionality was not changed by the GoogleTest commit;
- the new CPU test executable passed;
- the new HIP test executable has target-node evidence: 6/6 tests passed on Kunshan DCU;
- CUDA and KOKKOS remain unverified.

## 2. Git and remote state

| Item | State |
|---|---|
| Local branch | `master` |
| Local HEAD | `60681aaa` |
| Fork `origin/master` | `60681aaa` |
| `upstream/master` at handoff | `02899807` |
| Working tree | clean before adding this handoff document |
| Latest GoogleTest change on upstream | no PR opened |

The latest commit is pushed to the fork. For future upstream work, fetch
first, inspect the graph and changed files, then merge only after checking
ownership of `full.dat`, `solverstate`, residual outputs, and test infrastructure.

Recommended workflow:

```text
sync master
  -> create a feature branch
  -> make one logical change
  -> run the cheapest relevant tests
  -> commit and push to the fork
  -> open a PR to upstream master
  -> after merge, fetch and fast-forward local master
```

## 3. Residual regression architecture

Normal and strict residual regression use one architecture:

- one suite definition and runner;
- one canonical `test/baselines/residual-baseline.json`;
- one residual output path;
- one comparison implementation with mode-dependent tolerance and output precision.

The canonical database contains five cases and eight residual files. Normal mode
is intended for fast GitHub CI. Strict mode is intended for high-precision
validation on the Kunshan CPU queue. Strict mode uses an approximately
`1e-15` criterion and higher-precision residual output; it does not create
a second baseline system.

GoogleTest does not replace or duplicate this suite. It validates the lower-level
backend interface with small deterministic vectors. Full residual and MPI
regression remain integration-level checks.

## 4. Existing GPU/MPI evidence

The previous implementation has evidence for:

- CPU MPI: one-rank and 32-rank final hashes matched at tested sizes;
- HIP/DCU MPI: four ranks on four DCUs passed strict CPU/HIP comparison, with
  physical-state checks passing and four visible devices;
- one-card end-to-end speedups versus 32-rank CPU of approximately
  `3.85x`, `5.34x`, `6.80x`, and `9.70x`;
- four-card end-to-end speedups of approximately `0.57x`, `1.58x`,
  `5.39x`, and `13.10x`.

These are whole-process wall-clock comparisons. They do not establish CUDA or
KOKKOS compatibility.

## 5. GoogleTest change

The original `hello_test` remains as a GoogleTest smoke test. The new
contract test is
[tests/euler_backend_contract_test.cpp](../../../tests/euler_backend_contract_test.cpp).

It tests the backend directly through the existing interface:

- backend identity and accelerator metadata;
- `CreateState`, `Upload`, `Advance`, and `Download`;
- `FullTrace` and `NoTrace`;
- state reuse;
- finite and positive physical state;
- invalid problem, upload, step, and trace options;
- HIP device visibility when compiled with HIP.

The root CMake adds `oneflow_euler_backend_cpu_test`. The Kunshan port adds
`oneflow_1d_euler_hip_backend_test` only when both options are enabled:

```bash
-DONEFLOW_1D_ENABLE_GTEST=ON
-DONEFLOW_1D_ENABLE_HIP=ON
```

Both options are off by default in the Kunshan port. Developers without HIP or
DCU therefore do not need to build or run the GPU test.

Validated results:

- root-engine CPU backend contract: 5/5;
- standalone Kunshan-port CPU backend contract: 5/5;
- original `HelloTest`: 2/2.

The HIP target has now been compiled and run on a Kunshan DCU compute node.
The direct GoogleTest binary passed 6/6, and CTest now discovers and executes
all 6 `HIP.*` tests successfully after the backend-specific prefix and stale
discovery fixes.

## 6. Kunshan/DCU execution boundary

The `scnet-hpc` Kunshan profile identifies `kshdnormal` as the DCU
queue, requires at least `dcu:1`, and targets `gfx906`. Resource
requests must be generated from the profile and preserve its fixed
CPU/memory/GPU relationship. Do not invent a resource tuple.

HIP GoogleTest is a low-cost one-card compute-node test. Four cards remain for
MPI rank/device mapping and end-to-end integration regression.

A missing device or failed HIP runtime initialization must produce a nonzero
job status. A skipped test must not be reported as successful DCU validation.
Login-node package visibility is not compute-node compatibility evidence.

Earlier SCNet HTTP 503 and old-session API 400 failures were infrastructure or
session-context failures, not OneFLOW test results.

## 7. Next actions

1. Completed on 2026-09-03: submitted a profile-aligned one-card Kunshan DCU job.
2. Completed: built and ran `oneflow_1d_euler_hip_backend_test` on the compute-node
   HIP/DTK toolchain; direct GoogleTest passed 6/6.
3. Completed: recorded toolchain, architecture, visible-device evidence, exit
   status, and test summary in the maintained porting playbook.
4. Completed: added the reusable job script and sanitized validation notes.
5. Keep CPU GoogleTest in GitHub CI; run HIP GoogleTest only for HIP-related
   changes, scheduled validation, or an explicit release check.
6. Keep the `CPU.`/`HIP.` CTest prefixes and discovery-file refresh logic when
   extending the contract suite.
7. Split the contract test into lifecycle, numerical, and error-handling files
   only if it grows materially; the current test remains a coherent unit.

## 8. First files to inspect

- [backend contract test](../../../tests/euler_backend_contract_test.cpp)
- [root test CMake](../../../tests/CMakeLists.txt)
- [Kunshan port CMake](../../../ports/kunshan/oneflow_1d_hip/CMakeLists.txt)
- [Euler backend interface](../../../ports/kunshan/oneflow_1d_hip/OneDEulerBackend.h)
- [residual integration report](../regression/oneflow-residual-integration-20260901.html)
- [current performance report](../performance/oneflow-euler-performance-current.md)
- [report directory rules](../README.md)

## 9. Post-handoff DCU validation

The one-card HIP contract validation was completed on 2026-09-03 using the
profile-aligned `kshdnormal` request (1 node, 1 task, 8 CPU, `dcu:1`, 27G).
DTK 26.04 / Clang 17 compiled the HIP GoogleTest target for `gfx906`; the
compute-node binary passed 6/6 tests and the Slurm job completed with exit code
`0:0`. The reusable job is
[`ci/kunshan/euler-dcu-gtest.slurm`](../../../ci/kunshan/euler-dcu-gtest.slurm).

The run also captured two reusable rules: resolve CMake before `module purge`,
and never treat CTest's `No tests were found!` with exit code 0 as validation
success. The CTest harness is now fixed with backend-specific prefixes and a
fresh PRE_TEST discovery list; the broader migration guidance is maintained in
[`oneflow-dcu-porting-playbook-20260903.md`](oneflow-dcu-porting-playbook-20260903.md).
