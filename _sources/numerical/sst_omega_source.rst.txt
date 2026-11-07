SST omega-equation production, dissipation, and cross-diffusion
===============================================================

The specific-dissipation equation closes SST / Wilcox / EASM-kw families
with blended coefficients and an optional cross-diffusion term.

Blending
--------

.. math::

   \phi = F_1\phi_1 + (1-F_1)\phi_2

Menter production / dissipation
-------------------------------

.. math::

   P_\omega = \alpha\,\frac{\rho}{\mu_t} P_k Re,
   \qquad
   D_\omega = f_\beta\,\beta\,\rho\,\omega^2

.. math::

   CD_{k\omega}
   = 2(1-F_1)\rho\sigma_{\omega 2}
   (\nabla k\cdot\nabla\omega)/\omega

Wilcox-style production uses
:math:`P_\omega = \alpha\,(\omega/k) P_k` and
:math:`CD = \sigma_d\rho\max(\nabla k\cdot\nabla\omega,0)/\omega`.

References
----------

* Menter, F. R. Two-equation eddy-viscosity turbulence models (SST).
* Wilcox, D. C. Turbulence Modeling for CFD.
* NASA Turbulence Modeling Resource — SST
  (https://turbmodels.larc.nasa.gov/sst.html).
* See ``sst_k_source.rst`` and ``references.rst``.

Pure APIs: ``SstBlend``, ``SstProdwMenter``, ``SstProdwWilcox``,
``SstDissw``, ``SstCrossDiffMenter``, ``SstCrossDiffWilcox``.
Tests: ``turb_sst_policy_test``.
