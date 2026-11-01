# Contract Workbench knowledge catalog

Structured numerical knowledge for developers:

- taxonomy: family / genus (SA, SST, NS-...) / section
- plain formula + LaTeX source + optional mathtext preview
- principle, background, references
- Sphinx excerpt + open `.rst` / source files
- gated code sketch (checkbox required before copy)
- optional theory curve PNG from repo-relative CSV

All paths resolve from Workbench `repository_root`.

## Writable cache (do not use the source tree)

Rendered formula/curve previews are written under the user data root, never under
``tools/verification/data/``:

* default: ``~/.oneflow/contract-workbench/preview/``
* override root: environment variable ``ONEFLOW_HOME``

See ``tools/oneflow_gui/user_data_paths.py``.

## Plot style (axis labels)

Theory curve PNG generators under ``tools/verification/plot_*.py`` must:

* use descriptive **xlabel** / **ylabel** expressions with Matplotlib mathtext;
* put the governing **formula** in the title or subtitle;
* keep expressions within Matplotlib mathtext support; the preview is not a full LaTeX engine;
* call ``plot_style.apply_theory_style`` / ``save_figure``.

Do not commit bare-line plots without axes text. Runtime GUI previews still write under
user data paths (``~/.oneflow/contract-workbench/preview``), not the source tree.


## Generalized contract direction

This directory currently implements the numerical/formula specialization of the
Contract Workbench. A broader cross-domain model is under discussion in
[Verification Contract Workbench proposal](../../doc/plans/verification-contract-workbench.md).
Do not migrate the numerical catalog or add a generic runner until the shared
contract envelope, evidence semantics, and migration path have been reviewed.


## Read-only shared contract model (incremental migration)

`contract_catalog.py` is the first Qt-independent compatibility layer:

- `adapt_theory_entry(entry)` adapts the existing numerical `TheoryEntry`.
- `adapt_geometry_scenario(scenario)` adapts an already loaded Geometry
  scenario dictionary.
- `validate_contract_summaries(root, summaries)` checks stable IDs and
  repository-relative implementation, documentation, theory-data, and
  linked-test source paths.

The adapters do not modify either source catalog, execute tests, or persist
results. Evidence statuses such as `declared_not_executed` and
`linked_not_executed` are intentional: a reference to an API or GoogleTest
declaration is not proof that the test ran.

Run the headless compatibility tests from the repository root:

```powershell
python -m unittest tools.verification.test_contract_catalog
```

See [the staged architecture proposal](../../doc/plans/verification-contract-workbench.md).

## Combined read-only catalog

`contract_catalog.py` now exposes `load_contract_catalog(repository_root)`,
which adapts all existing Numerical entries and Geometry scenario JSON files
into namespaced immutable summaries:

- `numerical/<source-id>`
- `geometry/<source-id>`
- declared `runtime/<source-id>` and `case/<source-id>` contracts

The original ID is retained as `source_id`. The loader uses the existing
catalog/scenario readers, preserves the source formats, and returns diagnostics
for duplicate IDs, missing or escaping paths, stale linked GoogleTest
declarations, and unresolved Numerical API references (using the existing
Numerical validator). It does not execute tests or record results.

The Python test inventory discovers the new tests under
`tools/verification/test_contract_catalog.py`; this is covered by
`tools/oneflow_gui/test_test_inventory.py`.

Run focused checks:

```powershell
python -m unittest tools.verification.test_contract_catalog
python -m unittest tools.oneflow_gui.test_test_inventory
```

## Versioned cross-domain contract definitions

The first declared Runtime contract is stored at
`contracts/runtime/process_lifecycle.json`. It uses schema version 1 and
demonstrates that the common navigation envelope is not limited to numerical
formulas or geometry fixtures. Its payload describes a lifecycle state
transition, ownership boundaries, implementation paths, architecture, claims,
limitations, and typed GoogleTest evidence links.

`adapt_declared_contract` loads these versioned definitions into the same
immutable `ContractSummary` model. `load_contract_catalog` now combines:

- existing Numerical `TheoryEntry` adapters;
- existing Geometry scenario adapters;
- explicitly declared cross-domain JSON contracts.

A linked test remains `linked_not_executed`; catalog loading checks that its
GoogleTest declaration exists in the declared source file but does not run it
or infer a passing result. The runtime contract status `implemented` describes
implementation maturity, not verification completion. Its claims and
`does_not_prove` boundaries are part of the definition.

Run the headless catalog tests:

```powershell
python -m unittest tools.verification.test_contract_catalog
```

This is an initial runtime example, not a general runner or an invitation to
migrate all existing catalogs at once.

## First typed CFD case contract

The first Case-domain definition is
`contracts/case/laminar_flat_plate_slau2.json`. It describes the existing
`test/cfd/plateuns2dslau2` case selected by `test/suites/cpu-serial.txt`.
It records the solver, dimension, laminar Navier-Stokes model, SLAU2 flux,
50-step setup, input files, output comparison manifest, residual baseline,
and retained Blasius reference profiles.

The shared envelope remains common, while the case-specific fields are parsed
into the immutable `CaseContractData` payload. In addition to path and suite
checks, the catalog requires exact equality between the contract's output list
and `autotest/compare.txt`, checks that residual output channels exist in the
baseline, verifies the iteration columns and configured row count, and checks
the residual tolerance. The catalog validates that:

- the suite and referenced input/reference files exist;
- the case directory exists and is selected by the declared suite;
- every declared comparison output appears in the case's
  `autotest/compare.txt`;
- the residual database uses the supported schema, contains this case, and
  agrees with the declared absolute tolerance.

The analytical Blasius files are explicitly marked as references that the
current regression runner does not automatically compare. They must not be
mistaken for executed validation evidence. Likewise, the checked-in regression
baseline establishes reproducibility against an accepted result, not
independent physical validation or grid-convergence evidence.

This is the first domain-specific typed payload. It does not turn every case
into a contract: one case may provide evidence for several numerical,
geometry, boundary-condition, and runtime contracts. Nor does the catalog
launch `test.py`; execution remains with the established regression runner.


