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
CTest prefixes and the shared helper attaches
`hardware;hip;dcu` metadata. The cluster CTest 2.8 compatibility path exposes
`hardware` as the selectable label, so the DCU runner uses `-L hardware -R HIP`;
an empty discovery result or a zero exit code without the expected summary is
not successful validation.

Before configure, the runner executes [`dcu-preflight.sh`](dcu-preflight.sh).
It requires a successful `rocminfo` probe containing the requested
`ONEFLOW_HIP_ARCHITECTURES` agent (default `gfx906`) and a readable
`amd_comgr` CMake package directory. A node that cannot open `/dev/kfd`, has
the wrong architecture, or lacks `amd_comgr` now stops with
`preflight=1` in `result.txt`; it is not reported as a build or workload
failure.

The root project exposes the same HIP contract through
`-DONEFLOW_ENABLE_HIP_TESTS=ON`. It is deliberately opt-in: a normal CPU build
keeps the root CTest suite hardware-independent, while a target-node build must
also provide `CMAKE_HIP_COMPILER` and `CMAKE_HIP_ARCHITECTURES` (currently
`gfx906` on Kunshan Z100). The standalone project and the root project share the
same CMake registration helper.

## Root main-solver HIP/3D stage benchmark

For the production root solver (which remains CPU-only by default), use [`f3-main-solver-benchmark.slurm`](f3-main-solver-benchmark.slurm) on `kshdnormal` with `dcu:1`. The runner reuses a target-node root HIP build, checks the HIP smoke binary, non-empty root HIP GoogleTest, and hardware HIP CTest before running the 3D `m6wingroe_sa` stage gate. The trace gate compares legacy CPU, CPU batch, and HIP batch records and physical/semantic invariants before any timing is emitted. Set `ONEFLOW_F3_TRACE_STEPS` to choose the gate length; timing uses the same step count, one warmup, and three repeats by default. The measured wall-clock includes MPI launch, initialization, host packing, H2D, kernel, D2H, and output, so a result below `1.0x` is not an acceleration claim. The stream verifier is used for large traces to keep memory bounded.

For the current F2 vertical slice, [`f2-main-solver-stage.slurm`](f2-main-solver-stage.slurm) also supports `ONEFLOW_STAGE_BREAKDOWN_FILE` aggregation. The allocation may retain 8 CPU, 27G and `dcu:1`, but the checked-in single-zone `m6wingroe_sa` case must keep `ONEFLOW_F2_STAGE_BENCH_CPU_RANKS=1`. Set it to 8 only when the selected input is already partitioned into enough active zones. `ONEFLOW_F2_STAGE_TRACE_TIMEOUT_SECONDS` controls the short accuracy-trace timeout, while `ONEFLOW_F2_STAGE_BENCH_TIMEOUT_SECONDS` controls each warmup or formal repeat. The generated `breakdown.tsv` reads only `repeat-*` files and must contain no warmup rows.

Multi-rank main-solver timing is valid only for an already partitioned
multi-zone case. The stage, paired-timing, and stability runners verify that
every requested rank writes a breakdown file and records meaningful solver work.
A single-zone case run with multiple MPI processes therefore fails instead of
being reported as an N-rank CPU baseline. The checked-in `m6wingroe_sa` case is
single-zone, so its default CPU accuracy and stability basis is one rank.

Use [`f2-main-solver-stability.slurm`](f2-main-solver-stability.slurm) for a longer stability gate without writing the multi-gigabyte full-array trace. It prepares the same case for legacy CPU, CPU batch and HIP batch, runs them in that diagnostic order, and writes `stability.tsv` with ranks, workload exit code, diagnostic-hit count and pass/fail status. A scheduler `COMPLETED` state is not sufficient: inspect `stability.tsv`, each mode log and the workload exit code. The short 3-step trace remains the numerical and physical oracle; the longer runner checks whether the same configuration survives the requested step count without NaN/Inf/non-physical diagnostics.
By default the stability runner pins both `cflst` and `cfled` to `ONEFLOW_F2_STABILITY_CFL`, so the named CFL is a fixed value rather than only the ramp start. Set `ONEFLOW_F2_STABILITY_CFL_END` explicitly to request a ramp. Set `ONEFLOW_F2_STABILITY_LIMITER` to `-1`, `0`, `1`, or `2` to make the limiter choice explicit; the runner records all effective settings in `configuration.tsv`.
For an already partitioned multi-zone fixture, set `ONEFLOW_F2_STABILITY_PARTITIONED_CASE=1` and use the same rank count for `ONEFLOW_F2_STABILITY_CPU_RANKS` and `ONEFLOW_F2_STABILITY_HIP_RANKS`. The runner then forces `simutask=Solve`, selects `ONEFLOW_F2_STABILITY_GRID_FILE` (default `grid/m6wing-8.ofl`), and does not re-enter the partition workflow. An optional `ONEFLOW_F2_STABILITY_FIXTURE_METADATA` (or `ONEFLOW_MZ_FIXTURE_METADATA`) is validated before MPI starts and is recorded in `fixture-validation.txt`; leaving it unset records `MULTIZONE_FIXTURE_UNCHECKED`.

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

## Current Kunshan CPU/MPI toolchain (default)

All new Kunshan CPU, CPU-MPI, and root main-solver validation must use this
toolchain unless a test is explicitly marked as a historical reproduction:

| Component | Required selection | Module |
|---|---|---|
| Compiler | GCC 16.2.0 | `gcc/16.2.0` |
| MPI | OpenMPI 5.0.11 | `openmpi/5.0.11` |
| CMake | 4.4.3 | `cmake/4.4.3` |

Record `gcc --version`, `mpirun --version`, and `cmake --version` in each
run's environment evidence. For MPI runs, also record the resolved `mpirun`
and `libmpi` paths and the relevant `OMPI_MCA_*` variables. HCOLL must be
disabled for the standard MPI regression runs with
`OMPI_MCA_coll=^hcoll`; do not assume it is disabled merely because its
component is unavailable. Verify the effective environment in the run log.
The F2 stage and stability runners and the checked-in Actions configuration
template load this toolchain by default. The standalone DCU/HIP build runners
may retain a separately documented DTK-compatible host compiler/MPI stack;
that exception applies only to those HIP build/runtime workflows and must be
stated in the run evidence. Do not use that exception for CPU or CPU-MPI
results.

The 2026-09-13 results below are historical measurements from GCC 9.3.0 /
OpenMPI 4.1.5 / CMake 3.25.0. They remain useful as a record of that run but
are not the current environment specification and must not be used for new
validation claims.

Some older standalone Slurm runners still encode the former compiler and MPI
modules for reproducing dated measurements. Treat those as historical until
they are migrated; for new CPU/MPI runs use the current modules above and
capture the actual versions in `environment.txt`. In particular,
`f3-main-solver-benchmark.slurm` is pinned to its archived build/toolchain, and
`euler-dcu-gtest.slurm`, `euler-dcu-matrix.slurm`, and
`euler-dcu-mpi4-regression.slurm` are DTK/HIP-specific toolchain exceptions;
none of these establishes the CPU/MPI default.

## Historical CPU regression environment (kshcnormal, 2026-09-13)

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
- Kunshan standalone HIP contract (`dcu:1`, `gfx906`): 9/9 GoogleTest and 9/9 CTest passed on DTK 26.04.

## Verified DCU/HIP environment (kshdnormal)

Measured on 2026-09-13. The HIP targets live in the standalone
`ports/kunshan/oneflow_1d_hip` CMake project and are built separately from the
CPU solver with the DTK toolchain:

- module `compiler/dtk/26.04`, clang/clang++ 17
- `-DONEFLOW_1D_ENABLE_HIP=ON -DONEFLOW_1D_ENABLE_GTEST=ON`
- `-DCMAKE_HIP_ARCHITECTURES=gfx906` plus the DTK `hsa-runtime64_DIR` and
  `HSA_HEADER` cache entries used by `euler-dcu-gtest.slurm`

### Correctness

- HIP contract test: 9/9 passed (GoogleTest and CTest, `HIP.` prefix), including WENO5 CPU-oracle comparison.
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

## Standard workspace layout and test suites

The cluster workspace follows one root with fixed sub-directories so that
source, builds and run artifacts never mix. Concrete absolute paths stay in
the cluster-side `README.md`; the structure is:

```text
<workspace root>/
├── src/OneFLOW/              source archive for the revision under test
├── deps/metis-install/       self-built METIS 5.0.1
├── builds/<variant>/         solver-cpu | port-cpu | port-dcu | port-cpu-mpi | port-dcu-mpi
├── runs/<YYYYMMDD>/<suite>/  append-only run artifacts
├── work/                     run working directories
├── tmp/<jobid>/              job TMPDIR
├── probes/                   device/environment probes
└── archive/                  historical workspaces and failed builds
```

### Standard suites

| Suite | Partition | Resources | Content | Pass criteria |
|---|---|---|---|---|
| `cpu-regression` | `kshcnormal` | 16 CPU, 54G | five-case normal+strict; port CPU contract test | 5/5 normal, 5/5 strict, 5/5 contract |
| `dcu-single` | `kshdnormal` | 8 CPU, 27G, `dcu:1` | HIP contract test; stateful benchmark, four sizes | 9/9 contract; max abs error 0 and matching checksums for benchmark |
| `cpu-mpi` | `kshcnormal` | 32 ranks × 1 CPU | 32-rank CPU MPI benchmark, four sizes | exit 0, hashes valid |
| `dcu-mpi` | `kshdnormal` | 4 ranks × 8 CPU, `dcu:4` | 1-rank and 4-rank DCU MPI benchmark, four sizes | exit 0, hashes match `cpu-mpi`, `visible_devices=4` |

Common scale: `nx = 65536 / 262144 / 1048576 / 4194304`, `steps=100 repeats=2 warmup=1`.

The standalone 1D Euler MPI port is not a zone-based 3D solver. Its
`Partition(globalNx, rank, ranks, ...)` path repartitions the same global cell
array whenever the MPI rank count changes, and the halo ring uses those rank
neighbors. The rank-scaling matrix therefore keeps `global_nx`, initial
conditions, `steps`, `warmup`, and `repeats` fixed while testing CPU
`8/16/32 ranks` against HIP `1/2/4 ranks` (one rank per DCU). Each run must
report `local_ranks == requested ranks`, `decomposition_active_ranks ==
requested ranks`, a non-empty local partition, the expected halo-exchange
count, and a CPU/HIP final-hash match.

This rule does not apply to the root 3D solver. The checked-in `m6wingroe_sa`
case is single-zone, and the supported HIP capability requires
`ZoneState::nLocal == 1`; running that case with many MPI ranks would not
create valid per-rank zones. A valid 3D multi-rank experiment must first build
an already partitioned multi-zone fixture, then regenerate the
`ownerRanks`/interface exchange plan for the requested rank count and verify
that every rank owns a non-empty zone and writes solver work.

The external fixture can carry a small JSON topology manifest without adding
the partitioned mesh to Git. Set `ONEFLOW_MZ_FIXTURE_METADATA` when invoking
`f2-main-solver-multizone-trace.sh` or
`f2-main-solver-multizone-timing.sh`; once set, the runner requires
`ci/kunshan/verify-multizone-fixture.py` to accept the manifest before
launching MPI. The manifest schema is
`oneflow.multizone.fixture/v1` and records one non-empty zone owner per rank,
interface endpoints (`OwnedCell` to `GhostCell`), reciprocal neighbor pairs,
equation-major send/receive counts, and `HostStaged` or `DeviceDirect` memory
space. The run writes `fixture-validation.txt` alongside its timing/trace
artifacts. The private config may provide the same path as
`KUNSHAN_MZ_FIXTURE_METADATA`. Leaving both variables unset is allowed for legacy private fixtures,
but records `MULTIZONE_FIXTURE_UNCHECKED` and must not be described as an
ownership/topology validation. A copyable two-zone manifest template is kept
in [`multizone-fixture-metadata.example.json`](multizone-fixture-metadata.example.json);
the eight-zone fixture should expand the same fields to ranks `0..7`.

For a real prepared fixture, do not hand-edit endpoint IDs. Set
`ONEFLOW_MZ_FIXTURE_EXPORT_DIR` to a private, shared directory before the
solver's interface-update task runs. Each rank writes
`rank-<rank>.json` using the runtime `InterFace` and `FaceTopo` mappings
(including owned/ghost cell IDs and the configured equation count). After the
MPI run, merge those shards and run the normal verifier:

```bash
python3 ci/kunshan/merge-multizone-fixture.py \
  --shards "$ONEFLOW_MZ_FIXTURE_EXPORT_DIR" \
  --output "$ONEFLOW_MZ_FIXTURE_EXPORT_DIR/fixture-metadata.json"
python3 ci/kunshan/verify-multizone-fixture.py \
  --metadata "$ONEFLOW_MZ_FIXTURE_EXPORT_DIR/fixture-metadata.json" \
  --expected-ranks 8
```

The merger rejects missing ranks, duplicate zones, unpaired send/receive
endpoints, role/count mismatches, and any topology that would require
guessing from aggregate face counts. Set
`ONEFLOW_MZ_FIXTURE_EQUATIONS` when the solver is not five-equation Euler;
`ONEFLOW_MZ_FIXTURE_MEMORY_SPACE` defaults to `HostStaged` and must only be
set to `DeviceDirect` after a GPU-aware MPI probe has passed.

`lifecycle_*_ms` sums over `repeats`; only datasets recorded with the same
`repeats` value are comparable (see the erratum in
`doc/reports/performance/oneflow-euler-performance-20260913.md`).

### Run artifacts

Every run writes `environment.txt`, per-stage logs and `result.txt` under
`runs/<YYYYMMDD>/<suite>/`. Historical runs are never overwritten. Slurm
scripts and scheduler logs stay separate from the workspace root, under the
submitter's `scripts/` and `scripts/logs/` directories.

### Reusable job scripts

The four job scripts used for the 2026-09-13 runs are kept cluster-side. They
assume a single `ROOT` variable; when the workspace is reorganized, update
`ROOT` and the `builds/` and `runs/` sub-paths only — the module setup and
benchmark invocations stay unchanged.
