Viscous face spectral radius (overview)
=======================================

Purpose
-------

For Navier-Stokes time-step control, the inviscid face spectral radius

.. math::

   \lambda_{\mathrm{inv}}=\frac{1}{2}A\bigl(|v_n|+c\bigr)

must be supplemented by a **viscous** estimate \(\lambda_{\mathrm{vis}}\).
Viscous limits dominate in boundary layers and can force a much smaller
local \(\Delta t\). OneFLOW provides three pure policies in
``TimeStepPolicy``:

* **Model 1** — Blazek-style \(\max\) of momentum and energy diffusivities,
  scaled by \(A^{2}\) (see :doc:`viscous_spectral_radius_model1`).
* **Model 2** — length-scale form with cell/wall distance \(d\)
  (see :doc:`viscous_spectral_radius_model2`).
* **Model 3** — simplest \(A^{2}/(\mathrm{Re}\,\rho)\) placeholder
  (see :doc:`viscous_spectral_radius_model3`).

These are **engineering bounds**, not exact eigenvalues of the viscous
flux Jacobian.

Theory basis
------------

A diffusion mode \(u_t=\nu u_{xx}\) decays at rate \(\nu k^{2}\). On a
finite-volume mesh, \(k\) is estimated either from face/volume geometry
(\(k\sim A/V\), leading to \(A^{2}\) forms) or from an explicit length
\(d\) (leading to \(1/d\) forms). Non-dimensional NS equations place a
factor \(1/\mathrm{Re}\) on viscous fluxes, so \(\mathrm{Re}\) appears in
the spectral estimates.

Comparison
----------

.. list-table::
   :header-rows: 1
   :widths: 12 40 28 20

   * - Model
     - Core formula
     - Inputs
     - Role
   * - 1
     - \(A^{2}\cdot 2\max(c_1,c_2)/(\mathrm{Re}\rho)\)
     - \(\mu,\mu_t,\mathrm{Pr},\gamma,A\)
     - Full Blazek-type channel split
   * - 2
     - \(\mu_{\mathrm{tot}} A/(\rho\,d\,\mathrm{Re})\)
     - \(\mu_{\mathrm{tot}},d,A\)
     - Explicit length scale
   * - 3
     - \(A^{2}/(\mathrm{Re}\rho)\)
     - \(A,\mathrm{Re},\rho\)
     - Cheapest bound

Source map
----------

* ``ViscousFaceSpectralRadiusModel1`` / ``Model2`` / ``Model3``
* ``InviscidFaceSpectralRadius``, ``CellTimeStepFromSpectrum``,
  ``CombineInvAndVisTimeStep``
* Header: ``codes/ns/include/TimeStepPolicy.h``

References (primary and comparative)
------------------------------------

**Primary textbook (form of model 1)**

* J. Blazek, *Computational Fluid Dynamics: Principles and Applications*,
  Elsevier (2nd/3rd ed.). Chapter on temporal discretisation /
  determination of the maximum time step. Structured-grid viscous
  spectral radius (schematic)::

      Lambda_v ~ max( (4/3)(mu_L+mu_T)/rho ,
                      gamma (mu_L/Pr_L + mu_T/Pr_T)/rho ) * (Delta S)^2

  with a multiplier \(C\) (often 1--4) when combining with convective
  radii in \(\Delta t\). Unstructured variants sum face contributions
  with the same \(\max\) of momentum and energy diffusivities.
  ISBN / DOI vary by edition; see Elsevier monograph
  ``9780080445069`` (2nd) and later editions.

**Open-source codes (comparative, not OneFLOW proofs)**

* **SU2** — local time step from convective and viscous spectral radii;
  viscous contribution of the form \(\sum C\,(\mu/\rho)\,S^{2}\) along
  edges (see e.g. AIAA 2015-1949 and SU2 theory/solver documentation).
  Useful for terminology and scaling; coefficients differ.
* **Other FV RANS solvers** — similar \(\lambda_{\mathrm{vis}}\propto\mu S^{2}/\rho\)
  or \(A^{2}\) estimates appear widely; treat as implementation folklore
  unless tied to a specific published formula.

**Commercial manuals**

* Public ANSYS Fluent / STAR-CCM+ theory guides discuss CFL and time-step
  limits at a high level; they do **not** publish OneFLOW-equivalent
  internal spectral-radius coefficients. Use only for terminology, not
  as a correctness proof for OneFLOW formulas.

**What literature does and does not prove**

* Literature supports the *structure* of model 1 (momentum vs energy
  diffusivity, \(A^{2}\) or \((\Delta S)^{2}\) scaling, eddy viscosity).
* Exact numerical prefactors (\(2\), \(1/2\), choice of \(d\)) are
  **implementation choices** locked by OneFLOW unit tests, not unique
  theorems.
* Open/closed source manuals are for **traceability and comparison**,
  not permission to copy proprietary text or claim identical behavior.
