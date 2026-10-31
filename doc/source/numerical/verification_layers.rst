Verification layers for numerical contracts
===========================================

Purpose: keep models (especially turbulence) **self-evolving** without
relying only on full-case CFD-versus-experiment.

Layers (all required for a claim to "count")
--------------------------------------------

1. **Pure contract** — formula isolated in ``*Policy`` with ``@doc-*``.
2. **Unit test** — same expression locked by gtest (no shadow arithmetic).
3. **Intermediate check** — closed-form curve/table, published constant
   set, or independent reference implementation (e.g. SA
   :math:`f_{v1}(\\chi)`, SST :math:`F_1(\\mathrm{arg}_1)`).
4. **Documentation** — Sphinx page states the math, role in the pipeline,
   and short references (paper, NASA TMR, open manuals).
5. **Full case** (optional later) — mesh solve vs experiment/benchmark.

A comparison is valid only when **formula + data/reference + (optional)
plot** agree. Closed commercial manuals may be named; proprietary text
is not copied.

References policy
-----------------

See ``references.rst``. Prefer public anchors (NASA TMR, open guides,
classic papers). Extend the living list when a new pure API is extracted.
