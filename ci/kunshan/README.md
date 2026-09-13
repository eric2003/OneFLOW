# Kunshan regression workflow

The workflow in `.github/workflows/kunshan-regression.yml` is independent of
the existing Linux, Windows, and documentation workflows. It is intentionally
manual-only because it consumes Slurm resources and repository secrets are not
available to workflows from untrusted forks. It also stops before SSH setup
unless the selected workflow ref is the trusted `master`, `main`, or `hexin`
branch. This repository currently uses `master` as its default branch.

## Required GitHub Actions secrets

| Secret | Purpose |
|---|---|
| `KUNSHAN_SSH_HOST` | SSH host name |
| `KUNSHAN_SSH_PORT` | SSH port |
| `KUNSHAN_SSH_USER` | SSH user |
| `KUNSHAN_SSH_PRIVATE_KEY` | Private key accepted by the cluster |
| `KUNSHAN_SSH_KNOWN_HOSTS` | Pinned OpenSSH known-hosts line |
| `KUNSHAN_REMOTE_ROOT` | Dedicated remote CI root; runs are placed below `runs/` |
| `KUNSHAN_CI_CONFIG` | Absolute path to the private remote configuration file |

If any secret is absent, the private key is invalid, BatchMode SSH login
fails, Slurm commands are unavailable, or the remote configuration is
unreadable, the workflow stops before creating a run directory or submitting
a job and prints a specific Actions error.

## Remote configuration

Copy `config.example.sh` to a private path on the Kunshan shared filesystem,
fill in the target-cluster modules and dependency locations, and set
`KUNSHAN_CI_CONFIG` to that path.

The MPI fixture is kept outside Git because its partitioned grid is large.
`KUNSHAN_MPI4_CASE_ROOT` must point to a complete case directory containing:

- the four-zone grid;
- scripts configured to use that grid;
- an `autotest/` directory containing the accepted four-rank CPU baseline.

Serial and MPI baselines are deliberately separate. The serial job uses one canonical residual database and exercises both legacy and max-digits10 output modes. The current checked-in profile submits this suite to the CPU partition with `KUNSHAN_ACCEL_BACKEND=CPU`. Accelerator profiles reuse the same runner by selecting `HIP`, `CUDA`, or `KOKKOS` in their private configuration and must provide matching modules, partition/GRES requests, and target-node validation.

## One-card HIP contract validation

For the isolated Kunshan Euler backend contract, use
[`euler-dcu-gtest.slurm`](euler-dcu-gtest.slurm). It requests one DCU and
explicitly enables `ONEFLOW_1D_ENABLE_HIP` and `ONEFLOW_1D_ENABLE_GTEST`. The
job executes the GoogleTest binary and CTest, and records CMake, compiler, HIP
architecture, `rocminfo`, and Slurm evidence. CPU/HIP tests use `CPU.`/`HIP.`
CTest prefixes; an empty discovery result or a zero exit code without the expected
summary is not successful validation.

## Execution

Run **OneFLOW Kunshan Regression** from the Actions page and select:

- `cpu-serial`: the five cases in `test/suites/cpu-serial.txt`; the job runs both normal (`1e-8`) and strict (`1e-15`) residual profiles;
- `mpi4`: the private four-zone/four-rank M6 fixture;
- `all`: serial first, then MPI.

Each workflow run uses:

```text
KUNSHAN_REMOTE_ROOT/runs/<github-run-id>_<attempt>/
```

The source checkout, build, work directories, Slurm logs, accounting output,
and summary remain isolated under that directory. Logs and summaries are also
uploaded as a GitHub Actions artifact for 14 days. Each normal/strict run writes a machine-readable `backend-<mode>.txt` containing the selected backend, residual profile, tolerance, and exit status. This is the common reporting contract for CPU, HIP, CUDA, and Kokkos runs; an adapter is not considered validated until the report comes from its target compute node.

## Verified CPU regression environment (kshcnormal)

Measured on 2026-09-13, building the full solver (MPI + METIS + CGNS) and
running the five-case CPU serial suite in both residual profiles.

| Component | Selection | Notes |
|---|---|---|
| Compiler | GCC 9.3.0 | module `compiler/gcc/9.3.0` |
| MPI | OpenMPI 4.1.5 | module `mpi/openmpi/gcc-9.3.0/4.1.5` |
| CMake | 3.25.0 | module `compiler/cmake/3.25.0` |
| Python | 3.8.10 | module `python/3.8.10` |
| METIS | 5.0.1, built from source | tarball shipped with PHengLEI under `3rdparty/` |
| CGNS | 4.2.0 (prebuilt static) | `/public/software/mathlib/CGNS-4.2.0/src` |

Dependency wiring uses the variables the top-level CMake reads:
`MPI_HOME_INC`, `MPI_HOME_LIB`, `METIS_HOME_INC`, `METIS_HOME_LIB`,
`CGNS_HOME_INC`, `CGNS_HOME_LIB`.

### Pitfalls measured on this cluster

1. **GCC 7.3.1 (`compiler/devtoolset/7.3.1`) no longer builds the tree.** The
   solver uses `std::filesystem` and C++17 exception-specification deduction;
   both fail on GCC 7. Use GCC 9.3.0 or newer.
2. **The `compiler/gcc/11.2.0`, `12.2.0` and `13.3.0` modulefiles are broken on
   this cluster.** They prepend `/public/software/compiler/gcc/<ver>/bin`, which
   does not exist; the actual installations are at
   `/public/software/compiler/gcc-<ver>/bin`. Either fix the modulefiles or
   export the real path explicitly. `compiler/gcc/9.3.0` is correct.
3. **OpenMPI 4.x ships no C++ bindings.** `mpi.h` still pulls in the
   compatibility C++ headers while the library no longer implements them,
   producing `undefined reference to MPI::...` at link time. Configure with
   `-DCMAKE_CXX_FLAGS="-DOMPI_SKIP_MPICXX"`; the solver only uses the C API.
4. **METIS is not installed cluster-wide.** Build METIS 5.0.1 from the source
   tarball shipped with PHengLEI, using the same compiler as the solver.

### Measured results

- `normal` profile (`1e-8`): 5/5 cases passed, max absolute residual
  difference 4.9e-11.
- `strict` profile (`1e-15`, `ONEFLOW_RESIDUAL_TEST_OUTPUT=1`): 5/5 cases
  passed, max absolute difference 1.1e-17.
- Kunshan port CPU contract test (`-DONEFLOW_1D_ENABLE_GTEST=ON`): 5/5 passed.

## Verified DCU/HIP environment (kshdnormal)

Measured on 2026-09-13. The HIP targets live in the standalone
`ports/kunshan/oneflow_1d_hip` CMake project and are built separately from the
CPU solver with the DTK toolchain:

- module `compiler/dtk/26.04`, clang/clang++ 17
- `-DONEFLOW_1D_ENABLE_HIP=ON -DONEFLOW_1D_ENABLE_GTEST=ON`
- `-DCMAKE_HIP_ARCHITECTURES=gfx906` plus the DTK `hsa-runtime64_DIR` and
  `HSA_HEADER` cache entries used by `euler-dcu-gtest.slurm`

### Correctness

- HIP contract test: 6/6 passed (GoogleTest and CTest, `HIP.` prefix).
- Stateful benchmark, all four sizes: `final_max_abs_error = 0.000000` and
  identical CPU/HIP checksums.

### Performance (one DCU vs single-thread CPU advance, 100 steps, 2 repeats)

| nx | CPU ms/step | HIP ms/step | steady speedup | kernel ratio |
|---:|---:|---:|---:|---:|
| 65,536 | 11.05 | 0.071 | 155.9x | 99.3% |
| 262,144 | 46.77 | 0.227 | 206.2x | 99.7% |
| 1,048,576 | 192.68 | 0.836 | 230.4x | 99.9% |
| 4,194,304 | 881.44 | 3.274 | 269.2x | 100.0% |

The speedup column is a compute-path micro-benchmark relative to a single
serial CPU thread, not an application-level comparison. End-to-end numbers
against the 32-rank CPU MPI baseline remain in the performance reports.

### Node-level caveat

One allocated node (`e12r1n04`) returned `Unable to open /dev/kfd read-write:
Resource temporarily unavailable` and no HIP device, while a later allocation
on another node ran normally. A device probe (see `rocminfo` artefacts) is
therefore worth recording before blaming the build; if a node reports no
`gfx906` agent, exclude it and resubmit.
