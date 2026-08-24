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
│   ├── residual-baseline.json
│   └── residual-baseline-e15.json
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
python3 test.py \
  "mpirun -np 1" \
  "/path/to/OneFLOW" \
  suites/cpu-serial.txt \
  baselines/residual-baseline-e15.json
```

The fourth argument is optional for legacy suites. When supplied, the runner
compares every generated `res.full.dat` and `turbres.full.dat` row against the
versioned database. Enable those machine-regression files explicitly with
`ONEFLOW_RESIDUAL_TEST_OUTPUT=1`; normal solver runs leave them disabled and
continue writing the existing `res.dat`/`turbres.dat` files unchanged.
`ONEFLOW_RESIDUAL_DB` provides the same explicit database override.

The e-15 database covers all residual channels for every one of the 50 steps.
Iteration columns and row counts must match exactly; each residual value must
satisfy `abs(actual - reference) <= 1e-15`.

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
doc/reports/gpu-backend-delivery.md
```
