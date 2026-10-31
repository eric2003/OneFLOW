Turbulence face spectral radius
===============================

Implicit turbulence sweeps accumulate inviscid and diffusive spectral
contributions on each face.

Inviscid contribution
---------------------

.. math::

   \lambda_{\mathrm{inv}} = \tfrac12 A |v_n|

with mid-face relative normal velocity. Left/right matrix entries use
:math:`\tfrac12 A(\mp v_n^{L/R} - |v_n|)`.

SA effective viscosity
----------------------

.. math::

   \nu_{\mathrm{et}} = \mu_l/\rho + \tilde\nu

Diffusive contribution
----------------------

.. math::

   \lambda_{\mathrm{diff}}
   = \frac{1}{Re}\,\frac{1}{\sigma}\,\nu_{\mathrm{eff}}\,\frac{A^2}{V}

References
----------

* Spalart–Allmaras model spectral estimates for implicit integration.
* See ``turb_scalar_upwind.rst`` and ``spectral_radius.rst``.

Pure APIs: ``TurbSpectrumPolicy.h``. Tests: ``turb_spectrum_policy_test``.
