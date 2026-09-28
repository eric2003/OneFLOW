# OneFLOW tests

The `test/` directory contains test runners, suite definitions, test cases,
and checked-in numerical baselines. Generated reports and run artifacts do
not belong here.

## Directory policy

```text
test/
├── README.md
├── suites/
│   └── cpu-serial.txt
├── baselines/
│   ├── README.md
│   ├── build_residual_db.py
│   ├── verify_residual_db.py
│   └── residual-baseline.json
├── test.py
├── test.txt
└── <case directories>/
    ├── autotest/
    ├── grid/
    └── script/
```

- Add reusable suite lists under `test/suites/` with stable names.
- Keep accepted numerical references under each case's `autotest/` directory.
- Keep structured residual baselines under `test/baselines/`; regenerate them
  only when the accepted numerical result changes.
- Write generated logs, results, reports, and temporary files outside
  `test/`.
- Maintain project reports under `doc/reports/` using stable filenames.
- Do not create a new dated report for every regression run. CI run evidence
  belongs in workflow artifacts; durable conclusions update the stable report.

## CPU regression baseline

The canonical CPU baseline contains the five existing OneFLOW
cases:

| Case | Coverage | Reference outputs |
|---|---|---|
| `plateuns2dslau2` | 2-D laminar flat plate, SLAU2, viscous flux and wall boundary handling | `aero.dat`, `res.dat`, `flatplate_cf.dat`, `flatplateflow.dat`, `wallaero.dat` |
| `rae2822_roe_sa` | 2-D transonic airfoil, Roe flux and Spalart-Allmaras turbulence | `aero.dat`, `wallaero.dat`, `res.dat`, `turbres.dat` |
| `m6wingroe_sa` | 3-D wing, Roe flux and Spalart-Allmaras turbulence | `aero.dat`, `wallaero.dat`, `res.dat`, `turbres.dat` |
| `plateuns2dslau2_34950_35000` | continuation flat-plate case from the accepted restart | `aero.dat`, `res.dat`, `flatplate_cf.dat`, `flatplateflow.dat`, `wallaero.dat` |
| `turbplateuns2droe_sa` | turbulent flat plate, Roe flux and turbulence model | `aero.dat`, `wallaero.dat`, `res.dat`, `turbres.dat`, `turbplate_cf.dat`, `turbplateflow.dat` |

The checked-in files under each case's `autotest/` directory are the CPU
reference results. Accelerator implementations must first match these results
through the existing numeric comparison tolerance before performance results
are considered.

Run only this baseline suite from the `test` directory:

```bash
python3 test.py \
  "mpirun -np 1" \
  "/path/to/OneFLOW" \
  suites/cpu-serial.txt \
  baselines/residual-baseline.json \
  1e-15
```

The fifth argument is an optional residual absolute tolerance override. The
canonical database covers all residual channels for every one of the 50 steps.
Iteration columns and row counts must match exactly; strict residual validation
uses `abs(actual - reference) <= 1e-15`.
`ONEFLOW_RESIDUAL_DB` and `ONEFLOW_RESIDUAL_TOLERANCE` provide environment
variable equivalents. Set `ONEFLOW_RESIDUAL_TEST_OUTPUT=1` for the strict run;
the solver then writes max-digits10 values to the same `res.dat` and
`turbres.dat` paths used by the normal run.

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

The maintained delivery report is:

```text
doc/reports/architecture/oneflow-gpu-backend-delivery.md
```
