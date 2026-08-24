# Residual baseline databases

`residual-baseline-e15.json` is the active high-precision CPU residual database
for the `cpu-serial` suite. It contains every row and every residual channel
from `res.full.dat` and `turbres.full.dat` for all 50 solver steps.

The solver produces these machine-regression files only when explicitly enabled:

```bash
export ONEFLOW_RESIDUAL_TEST_OUTPUT=1
```

Normal runs do not create `.full.dat` files. Existing `res.dat` and
`turbres.dat` formatting and behavior remain unchanged.

The high-precision comparison requires identical iteration columns and row
counts, rejects NaN/Inf, and applies an absolute tolerance of `1e-15` with no
relative tolerance. Values are written with `max_digits10` (17 significant
digits for OneFLOW's `double` `Real`).

Each database entry records the source reference and SHA256, variable names,
the complete residual curve, summaries, comparison policy, suite, launcher,
platform, and validated revision metadata.

Regenerate the active database only after accepting new CPU reference results:

```bash
python3 test/baselines/build_residual_db.py \
  --high-precision \
  --validated-commit <validated-revision> \
  --output test/baselines/residual-baseline-e15.json
```

Validate database structure, source hashes, and stored curves:

```bash
python3 test/baselines/verify_residual_db.py \
  test/baselines/residual-baseline-e15.json
```

`residual-baseline.json` is retained as the legacy low-precision database. It
uses the historical `res.dat`/`turbres.dat` files and is not the e-15 gate.
Per-run logs and scheduler records remain CI artifacts and do not create new
baseline files.
