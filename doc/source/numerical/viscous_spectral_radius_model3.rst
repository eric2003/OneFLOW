Viscous face spectral radius (model 3)
======================================

Purpose
-------

Model 3 is the cheapest viscous spectral placeholder:

.. math::

   \lambda_{\mathrm{vis}}=\frac{A^{2}}{\mathrm{Re}\,\rho}

Family overview: :doc:`viscous_spectral_radius`.

Formula and rationale
---------------------

Only geometry (:math:`A`), density, and Reynolds number enter. Dynamic
viscosity is **not** explicit: viscous magnitude is assumed absorbed into
the non-dimensional \(\mathrm{Re}\) scaling. This is useful as a robust
fallback or when a full \(\mu,\mathrm{Pr}\) field is unavailable.

It still scales as :math:`A^{2}\), consistent with a \(k\sim A/V\) diffusion
argument, but without the \(\max(c_1,c_2)\) structure of model 1.

Source map
----------

* ``ViscousFaceSpectralRadiusModel3`` in ``TimeStepPolicy``
* Existing unit test:
  ``TimeStepPolicy.ViscousFaceSpectralRadiusModel3IsArea2OverReRho``
* Contract Workbench: add when an intermediate curve entry is introduced

What this proves / does not prove
---------------------------------

**Proves:** \(A^{2}/(\mathrm{Re}\rho)\) algebra matches the pure API.

**Does not prove:** adequacy for high-\(\mu_t\) boundary layers (model 1
is preferred when viscosities are known).

References
----------

* Degenerate limit of area-squared viscous estimates when the
  diffusivity prefactor is absorbed into \(\mathrm{Re}\) (engineering
  simplification; not a separate classical theorem).
* Blazek (same time-step chapters) motivates *including* viscous spectral
  radii; model 3 is OneFLOW's minimal realization of that requirement.
* Family overview: :doc:`viscous_spectral_radius`.
