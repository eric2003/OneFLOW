# Residual baseline database

`residual-baseline.json` is the versioned CPU residual database for the
`cpu-serial` suite. It contains the complete `res.dat` and `turbres.dat`
curves, not only the final residual value.

Each entry records:

- source reference file and SHA256;
- variable names and every `(iter, sub-iter, residual...)` row;
- first/last values and per-channel maximum magnitude;
- comparison tolerances and index-column policy;
- suite, launcher, platform and validated commit metadata.

Regenerate only when an accepted numerical baseline changes:

```bash
python3 test/baselines/build_residual_db.py \
  --validated-commit <commit-that-passed-the-baseline>
```

Validate database structure, source hashes and stored curves:

```bash
python3 test/baselines/verify_residual_db.py
```

Normal CI runs read the database and compare generated residual curves against
it. Per-run logs and scheduler records remain CI artifacts; they do not create
new database files.
