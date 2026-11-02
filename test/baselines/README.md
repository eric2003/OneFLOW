# Residual baseline database

`residual-baseline.json` is the canonical CPU residual database for the
`cpu-serial` suite. It covers all five serial CPU cases and stores the
full-precision `res.dat`/`turbres.dat` curves. The same reference supports
normal (`1e-8`) and strict (`1e-15`) comparison profiles.

The solver keeps the historical five-digit output by default. CI can opt into
canonical full-precision serialization with:

```bash
export ONEFLOW_RESIDUAL_TEST_OUTPUT=1
```

The flag changes only residual text serialization; it does not change solver
arithmetic or backend selection. Full-precision mode writes the same
`res.dat`/`turbres.dat` paths, so ordinary and strict tests use one output
interface.

The strict comparison requires identical iteration columns and row counts,
rejects NaN/Inf, and applies an absolute tolerance of `1e-15` with no relative
tolerance. Normal comparison overrides the same database to `1e-8`.

Each database entry records the source reference and SHA256, variable names,
the complete residual curve, summaries, comparison policy, suite, launcher,
platform, and validated revision metadata.

Regenerate the canonical database only after accepting new CPU reference results:

```bash
python3 test/baselines/build_residual_db.py \
  --suite test/suites/cpu-serial.txt \
  --validated-commit <validated-revision> \
  --absolute-tolerance 1e-15 \
  --output test/baselines/residual-baseline.json
```

Validate database structure, source hashes, and stored curves:

```bash
python3 test/baselines/verify_residual_db.py \
  test/baselines/residual-baseline.json
```

Per-run logs and scheduler records remain CI artifacts and do not create new
baseline files.
