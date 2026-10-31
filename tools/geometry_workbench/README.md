# Geometry contract workbench

Table-level scenarios for interface dual-map / split contracts. Shared by:

* CLI: `python tools/geometry_workbench/__main__.py`
* OneFLOW Workbench GUI: home -> **Geometry Contracts...**

Not a GoogleTest launcher. The current stub is fixture-backed: it presents a
known table result without calling the production geometry algorithm. In the
GUI, both input tables and expected tables are editable JSON, with separate
restore buttons for the scenario defaults. Running validates the input JSON,
compares the fixed mock result against the edited expected contract, and renders
the mismatch and topology view. The current stub does not calculate output from
edited input; a regression test makes this boundary explicit. This makes the
scenario-editing, comparison, and visualization workflow testable before a
native backend exists. Editing values tests the workbench workflow and contract
checker, not the CFD algorithm. The `registered_child_partial_match` scenario
links to the authoritative C++ regression
`FaceSearchBasicTest.RegisteredChildResolutionPreservesCandidateIndexForPartialMatch`,
but the GUI does not yet execute that test or pass edited parameters to C++.
Native execution must be added as a separate adapter rather than duplicating the
resolver algorithm in Python.

## Optional static graphics

```bash
python __main__.py --plot-file /tmp/dual_map.png
```

Requires matplotlib. Uses the Agg backend (no GUI window). Skip this flag for
text-only runs in CI.

## Domains beyond geometry

Scenario JSON may set `"domain": "geometry" | "turbulence" | "solver" | ...`.
The workbench UI is shared; only the runner differs:

* geometry fixtures (current stub / future native C++)
* `analytic_*` intermediate formula tables (see `analytic_contracts.py`)
* future: launch matching GoogleTest with parameters

Example intermediate scenario: `scenarios/sa_fv1_chi_curve_intermediate.json`.
