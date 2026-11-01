# Intermediate numerical verification (Contract Workbench)

## Role

Developer-facing knowledge browser for numerical formulas used in OneFLOW:

- taxonomy: family / model (SA, SST, ...) / section (closure, source, ...)
- formula, principle, background, references
- optional theory CSV + static plot
- links to pure APIs and Sphinx pages

Paths are **repository-relative** (`tools/verification/data/...`), resolved from
the Workbench `repository_root`. Do not hard-code machine-specific absolute paths.

## GUI

`oneflow_gui` home: **Verification Curves...**
(window title: Contract Workbench — Formula & Theory Curves)

## Long-term

Folder name may become `tools/contract_workbench`. Catalog entries are a step
toward formula-driven documentation and, later, assisted code generation from
the same structured knowledge — still requiring human review and tests.
