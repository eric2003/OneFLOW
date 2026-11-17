Viscous face spectral radius (model 2)
======================================

Purpose
-------

Model 2 estimates \(\lambda_{\mathrm{vis}}\) using an explicit length
scale :math:`d` (cell size or wall distance). Family overview:
:doc:`viscous_spectral_radius`.

Formula
-------

.. math::

   c_1=\frac{2\mu_{\mathrm{tot}}}{\rho\,d\,\mathrm{Re}},\qquad
   \lambda_{\mathrm{vis}}=\frac{1}{2} c_1 A
   =\frac{\mu_{\mathrm{tot}} A}{\rho\,d\,\mathrm{Re}}

Derivation (sketch)
-------------------

Diffusion rate :math:`\nu/d^{2}` with
:math:`\nu\sim\mu/(\rho\,\mathrm{Re})` yields an inverse-time scale.
Mapping onto a face contribution linear in :math:`A` produces the
form above (engineering factors :math:`2` and :math:`1/2`).

Unlike model 1, there is no separate energy Prandtl channel; a single
total viscosity is used.

Source map
----------

* ``ViscousFaceSpectralRadiusModel2`` in ``TimeStepPolicy``
* Tests: ``tests/solver/timestep_policy_test.cpp`` (when intermediate
  curve tests are present)
* Contract Workbench id: ``viscous_spectral_radius_model2`` (when wired)

What this proves / does not prove
---------------------------------

**Proves (when tests present):** algebra and :math:`1/d` scaling versus
the pure API.

**Does not prove:** optimal geometric definition of :math:`d` on a mesh.

References
----------

* Diffusion time-scale argument (standard CFD notes; Blazek discusses
  viscous spectral limits in boundary layers in the same chapters as
  model 1).
* Length-based viscous CFL estimates appear in many unstructured codes;
  SU2-style \(\mu S^{2}/\rho\) is a related but **not identical** form
  (area-squared without explicit \(d\)).
* Family overview: :doc:`viscous_spectral_radius`.
