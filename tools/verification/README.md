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
