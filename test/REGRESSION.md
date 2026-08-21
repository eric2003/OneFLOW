# CPU regression baseline

The initial accelerator-development baseline contains three existing OneFLOW
cases:

| Case | Coverage | Reference outputs |
|---|---|---|
| `plateuns2dslau2` | 2-D laminar flat plate, SLAU2, viscous flux and wall boundary handling | `aero.dat`, `res.dat`, `flatplate_cf.dat`, `flatplateflow.dat`, `wallaero.dat` |
| `rae2822_roe_sa` | 2-D transonic airfoil, Roe flux and Spalart-Allmaras turbulence | `aero.dat`, `wallaero.dat`, `res.dat`, `turbres.dat` |
| `m6wingroe_sa` | 3-D wing, Roe flux and Spalart-Allmaras turbulence | `aero.dat`, `wallaero.dat`, `res.dat`, `turbres.dat` |

The checked-in files under each case's `autotest/` directory are the CPU
reference results. Accelerator implementations must first match these results
through the existing numeric comparison tolerance before performance results
are considered.

Run only this baseline suite from the `test` directory:

```bash
python3 test.py "mpirun -np 1" "/path/to/OneFLOW" regression_cpu.txt
```

The repository also contains a manual-only Kunshan workflow:

```text
.github/workflows/kunshan-regression.yml
```

It performs an SSH/Slurm authorization preflight before creating a remote
directory or submitting a job. If the required Kunshan secrets, pinned host
key, private key, or remote configuration are unavailable, the workflow stops
with an explicit error. See [`ci/kunshan/README.md`](../ci/kunshan/README.md)
for the required secrets and private cluster configuration.

The backend framework defaults to CPU and does not redirect numerical
operators yet. Runtime selection can be requested with
`ONEFLOW_ACCEL_BACKEND`; requesting an adapter that has not been built fails
with a clear error rather than silently falling back to CPU. The CMake cache
also accepts `HIP`, `CUDA`, and `KOKKOS` as reserved backend names so future
adapters can be added without changing the solver-facing interface.
