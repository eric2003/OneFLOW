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
checker, not the CFD algorithm. Native backend integration is reserved for
`RunSymmetricInterfaceSplitMappingPipeline`.

## Optional static graphics

```bash
python __main__.py --plot-file /tmp/dual_map.png
```

Requires matplotlib. Uses the Agg backend (no GUI window). Skip this flag for
text-only runs in CI.
