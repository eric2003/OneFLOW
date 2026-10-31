Spalart-Allmaras source residual
================================

After closure functions are known, the SA residual combines production,
destruction, and a limited cross-diffusion term.

.. math::

   P = c_{b1}(1-f_{t2})\tilde S\,\tilde\nu

.. math::

   D = (c_{w1} f_w - c_{b1} f_{t2})\,\tilde\nu\, r_s

.. math::

   \mathrm{Diff}
   = \min\big(r_p P,\,
   (c_{b2}/\sigma)|\nabla\tilde\nu|^2/Re\big)

.. math::

   S = P - D + \mathrm{Diff}

Spectral contribution uses the negative part of linearized
production-destruction rates:
:math:`-\tfrac12(x-|x|)`.

References
----------

* Spalart, P. R. and Allmaras, S. R. A one-equation turbulence model.
* NASA Turbulence Modeling Resource (SA).
* See ``sa_closure.rst`` and ``references.rst``.

Pure APIs: ``SaProduction``, ``SaDestruction``, ``SaDiffusion``,
``SaLimitedDiffusion``, ``SaSourceResidual``, ``SaNegativePart``.
Tests: ``turb_sa_policy_test``.
