# AI task template (OneFLOW hexin geometry / architecture)

Use this shape instead of long free-form Chinese/English paragraphs.
Machines parse **fields**; humans still read **Goals**.

```text
REPO: OneFLOW
BRANCH: hexin
BASE: latest origin/hexin (or paste commit SHA)
LOCAL_TESTS: N/N passed

GOAL:
  One sentence outcome (e.g. "lock pure prefix of line interface split").

CONSTRAINTS:
  - tests only under tests/geometry
  - do not edit codes/**/CMakeLists.txt
  - source ASCII only
  - no silent product-behavior change without tests
  - prefer closing an architecture-map row over renaming stable helpers

ROI_PRIORITY (pick one primary):
  [ ] 1 solver/numerical contracts (limiter, Marray, residual path)
  [ ] 2 geometry HARD EDGE (locator + nodeMesh, not more pure predicates)
  [ ] 3 workbench native bind to table pipeline
  [ ] 4 full CGNS e2e
  [ ] 5 docs/diagrams/tags only

DO:
  - code change list (files)
  - tests added/updated
  - doc/source update (algorithm + diagram/math if applicable)
  - optional: @doc-* tags on pure header APIs

DO_NOT:
  - average effort across low-ROI pure helpers
  - strip orphan inverse rows without migration design
  - second algorithm implementation in Python/GUI

DELIVER:
  - zip with full relative paths of changed sources
  - BASE_COMMIT.txt + COMMIT_MESSAGE.txt
  - expected test count delta
```

## Diagrams and math in docs

- Flow: `.. code-block:: text` ASCII or future graphviz; label pure vs hard-edge.
- Formulas: Sphinx `.. math::` (mathjax enabled in `doc/source/conf.py`).
- Icons: use words in tables (Closed / Still open); avoid emoji in C++ sources.

## Comments -> docs

Stable tags in headers (`@doc-algorithm`, `@doc-stage`, `@doc-pure`, `@doc-formula`)
are the extractable surface. Narrative stays in `doc/source/`.
