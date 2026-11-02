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

Pure APIs: ``SaEddyViscosity`` and shared ``TurbCappedEddyViscosity`` (see ``turb_viscosity_policy_test``).
Wiring: ``TurbSrcFlux::CalcCellVist1Equ``.


Density validation in the SA field loop is performed on the current cell value
read from the mean-flow density field. The validity predicate is covered
independently by ``turb_viscosity_policy_test``. This ordering matters: checking
the shared ``turbcom.rho`` before loading the cell value would inspect stale
state (or an uninitialized value) rather than the cell currently being processed.


Complete per-cell closure
-------------------------

``SaComputeEddyViscosity(input)`` now composes the complete SA eddy-viscosity
calculation from primitive cell state: density, working variable, laminar
viscosity, :math:`c_{v1}^3`, and the configured viscosity-ratio limit. Its
result exposes :math:`\chi`, :math:`f_{v1}`, the uncapped eddy viscosity,
and the final capped value. The production adapter
``TurbSrcFlux::CalcCellVist1Equ`` maps the solver state into this input and
retains the intermediate :math:`\chi`, :math:`\chi^3`, and
:math:`f_{v1}` values used by the surrounding solver state.

The focused tests verify both the full intermediate chain and the shared
laminar-viscosity cap. This is a cell-closure contract; it does not replace
the field-loop and mesh-level verification described above.
