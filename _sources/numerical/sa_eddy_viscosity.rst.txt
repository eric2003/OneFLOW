SA eddy viscosity and intermediate verification
==============================================

The turbulent eddy viscosity for Spalart-Allmaras is recovered from the
working variable through :math:`f_{v1}`:

.. math::

   \mu_t = \rho\,\tilde\nu\,f_{v1}(\chi),
   \qquad
   \chi = \tilde\nu/\nu

with an optional solver cap
:math:`\mu_t \leftarrow \min(\mu_t,\, r_{\max}\mu_l)`.

Intermediate verification (not a full-run comparison)
-----------------------------------------------------

Unit tests alone prove the algebra. For turbulence models it is also
useful to lock **intermediate curves** against the closed-form theory
(and, later, against published calibration tables or independent codes):

* **Theory**: evaluate :math:`f_{v1}(\chi)` on a fixed :math:`\chi`
  grid and compare to the analytic expression (see
  ``TurbSaPolicy.Fv1ChiCurveMatchesClosedForm``).
* **Constants**: use documented SA values (:math:`c_{v1}=7.1`, NASA TMR).
* **Future**: optional CSV/plot of :math:`f_{v1}`, :math:`f_w(r)`, wall
  law profiles against Spalart–Allmaras / TMR reference data — still
  driven by the same pure APIs, not by a full mesh solve.

This sits between pure unit tests and end-to-end CFD-vs-experiment
comparisons: formulas + references + (optional) plots must agree.

References
----------

* Spalart, P. R. and Allmaras, S. R. A one-equation turbulence model.
* NASA Turbulence Modeling Resource — Spalart–Allmaras
  (https://turbmodels.larc.nasa.gov/spalart.html).
* See ``sa_closure.rst``, ``sa_source.rst``, and ``references.rst``.

Pure APIs: ``SaEddyViscosity``, ``SaCappedEddyViscosity``.
Wiring: ``TurbSrcFlux::CalcCellVist1Equ``.
