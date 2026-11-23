# OneFLOW verification

Placeholder until the audit selects a real numerical/regression case.

`run_all.sh` exits 2 (`UNVERIFIED`), so an integration gate refuses it.
Use OneFLOW's existing `test/` regression suite and `tests/` GoogleTest contracts
as the source of truth; do not invent a new numerical baseline.
