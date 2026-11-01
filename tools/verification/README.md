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

* set **xlabel** / **ylabel** with plain ASCII names (readable without math fonts);
* put the **formula** in the title subtitle;
* call ``plot_style.apply_theory_style`` / ``save_figure``.

Do not commit bare-line plots without axes text. Runtime GUI previews still write under
user data paths (``~/.oneflow/contract-workbench/preview``), not the source tree.
