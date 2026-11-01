Turbulence face spectral radius
===============================

Inviscid contribution
---------------------

.. math::

   \lambda_{\mathrm{inv}} = \frac{A}{2}|v_n|

Pure API: ``TurbInviscidFaceSpectralRadius``. Intermediate table:
``tools/verification/data/turb_inviscid_spectrum_theory.csv`` (reference face area
:math:`A=2`).

.. image:: figures/turb_inviscid_spectrum_theory.png
   :alt: Turb inviscid spectrum vs abs vn
   :width: 480

Diffusive contribution (SA)
---------------------------

.. math::

   \nu_{\mathrm{et}} = \frac{\mu_l}{\rho}+\tilde\nu,\qquad
   \lambda_{\mathrm{diff}} = \frac{1}{Re}\frac{1}{\sigma}\nu_{\mathrm{et}}
   \frac{A^{2}}{V}

Locks: ``turb_spectrum_policy_test``. Contract Workbench entry
``turb_inviscid_spectrum``.

**Proves:** pure spectral-radius algebra. **Does not prove:** turb LUSGS case
convergence.
