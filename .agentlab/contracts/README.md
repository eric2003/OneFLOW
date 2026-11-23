# OneFLOW contracts

The generic CFD contracts in this directory are **starter templates**, not OneFLOW contracts.
They describe a small `ISpatialOperator`/`ITimeIntegrator` boundary that OneFLOW does not
currently expose in that form.

During the audit:

- treat `interfaces.md`, `data_layout.md`, and `core_api.md` as examples only;
- derive real boundaries from `codes/main`, `codes/io`, `codes/ns`, `codes/residual`,
  `codes/grad`, `codes/solver`, `codes/time`, `codes/unsteady`, `codes/update`, and their tests;
- propose replacement contracts in the audit report or an RFC.

Do not implement against the generic contracts without a human-approved replacement.
