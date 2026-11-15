Verification layers (copyable policy)
=====================================

Goal
----
Evolve numerical models with contracts that remain testable without
depending only on full-case CFD-versus-experiment.

Layers
------
1. Pure contract — formula in ``*Policy`` (``@doc-*`` tags).
2. Unit test — gtest locks the same expression.
3. Intermediate check — closed-form curve/table, published constants,
   or independent reference code (examples: SA \(f_{v1}(\chi)\),
   SST \(F_1(\mathrm{arg}_1)\), Wilcox \(f_w(x_w)\)).
4. Documentation — math, pipeline role, short bibliography.
5. Full case (optional) — mesh solution versus experiment/benchmark.

Acceptance rule
---------------
A claim is valid only when **formula + reference data/source +
(optional) plot** agree. Proprietary manuals may be named; proprietary
text is not copied.

References
----------
Prefer public anchors (papers, NASA TMR, open guides). See
``references.rst``. Extend the list when a new pure API is extracted.
