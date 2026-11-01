# Intermediate numerical verification

## Product name: Contract Workbench

Shared operator surface inside `oneflow_gui` for:

| Domain | Tooling today | GUI entry |
|--------|---------------|-----------|
| Geometry table contracts | `tools/geometry_workbench` | Geometry Contracts... |
| Turbulence theory curves | `tools/verification` | Verification Curves... |

Do **not** rename packages in this step. Long-term folder name candidate:
`tools/contract_workbench` (geometry + numerical domains).

Theory curves are closed-form tables/plots. They are not production C++ and
not a substitute for GoogleTest or case regression.

# Intermediate verification

See sa_closure_curves.py and plot_sa_fw.py.
| `plot_sa_fv2.py` | Optional PNG of fv2(chi) |
