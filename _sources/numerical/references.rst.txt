References and external manuals
===============================

OneFLOW numerical docs cite sources when a formula is standard CFD theory
or when an external manual clarifies design choices. Citations are **not**
required for every algebraic rewrite of an already-locked pure API.

When to cite
------------

* Classic closed-form laws (Sutherland, ideal-gas relations, Roe averages).
* Named schemes (Roe, AUSM, Harten entropy fix) and parameter defaults.
* Spectral-radius / CFL estimates that match textbook presentations.
* Optional comparison notes to open manuals (OpenFOAM, SU2, Code_Saturne)
  when explaining *why* a pure contract is shaped a certain way.

How to cite
-----------

* Prefer a short **References** subsection at the end of the relevant
  ``doc/source/numerical/*.rst`` page (author, title, edition/year).
* Public URLs are welcome for open guides; avoid pasting long copyrighted
  excerpts.
* Closed commercial theory manuals may be named at a high level when the
  team already relies on them; do not reproduce proprietary text.
* Code ``@doc-*`` tags stay formula-local; bibliography lives in Sphinx.

Living list (seed)
------------------

* Blazek, J. *Computational Fluid Dynamics: Principles and Applications*.
* Toro, E. F. *Riemann Solvers and Numerical Methods for Fluid Dynamics*.
* Roe, P. L. Approximate Riemann solvers (J. Comput. Phys.).
* Liou, M.-S. AUSM-family schemes.
* OpenFOAM User Guide / Programmer's Guide (open).
* SU2 documentation (open).

This page is a seed; extend it when a pure-API doc needs a durable citation.

Self-evolving documentation loop
--------------------------------

The numerical subsystem is meant to evolve as a closed loop:

1. **Pure contract** in ``*Policy.h`` (``@doc-*`` tags, formulas).
2. **gtest** locks the same expressions (no shadow arithmetic).
3. **Sphinx page** states the algorithm, pipeline role, and short refs.
4. **Wiring** replaces hand-rolled copies in solvers (flux, LUSGS, CFL).

External manuals (open or internal) enter at step 3 as *anchors*, not as
paste. When a new pure API is extracted, prefer adding one citation if the
formula is named in the literature; otherwise link the sibling policy page.

* Spalart, P. R. and Allmaras, S. R. A one-equation turbulence model for
  aerodynamic flows.
* NASA Turbulence Modeling Resource (https://turbmodels.larc.nasa.gov/).
