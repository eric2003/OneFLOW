# Geometry contract workbench

Table-level scenarios for interface dual-map / split contracts. Shared by:

* CLI: `python tools/geometry_workbench/__main__.py`
* OneFLOW Workbench GUI: home → **Geometry Contracts...**

Not a GoogleTest launcher. Stub backend returns expected tables; native backend
is reserved for `RunSymmetricInterfaceSplitMappingPipeline`.

## Optional static graphics

```bash
python __main__.py --plot-file /tmp/dual_map.png
```

Requires matplotlib. Uses the Agg backend (no GUI window). Skip this flag for
text-only runs in CI.
