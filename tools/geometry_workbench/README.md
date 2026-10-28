# Geometry workbench (table scenarios)

Development aid for interface dual-map / split contracts. **Not** a GoogleTest
runner and **not** a second implementation of the C++ algorithms.

## Matplotlib vs interactive backends

| Backend | Role | Interactivity |
|---------|------|----------------|
| ASCII | Always on | None (copy/paste friendly) |
| matplotlib | Optional sketch | Toolbar zoom/pan; poor for editing |
| plotly | Optional sketch | Browser zoom/hover/pan; better exploration |

**Editable inputs and pass/fail judgment live in the host (CLI JSON today,
future panel), not inside the plot library.** Matplotlib is fine for a quick
picture; it is not the product UI. Prefer plotly when exploring graphs; use a
real form/grid for changing `parents` / `children` / `l2g`.

## Run

```bash
cd tools/geometry_workbench
python __main__.py
python __main__.py --plot none
python __main__.py --plot plotly    # if plotly installed
```

Optional: `pip install plotly` or `matplotlib`. Core path needs only stdlib.

## Backend swap

- `backend=stub`: returns `expected` from the JSON (UI track can ship now).
- `backend=native`: reserved for binding `RunSymmetricInterfaceSplitMappingPipeline`
  (C++). Do not re-implement split logic in Python.
