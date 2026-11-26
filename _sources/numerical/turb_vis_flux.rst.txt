Turbulence scalar diffusive face flux
=====================================

SA and two-equation models share the same scalar diffusion write-back:

.. math::

   F = -\frac{1}{Re}\,\Gamma\,\frac{\partial\phi}{\partial n}\,A

Effective diffusivity :math:`\Gamma`:

* SA: :math:`(\nu_l+\tilde\nu)/\sigma`
* SST: :math:`\mu_l+\sigma_k\mu_t` (and :math:`\sigma_\omega` for omega)

References
----------

* Spalart–Allmaras; Menter SST (diffusion coefficients).
* See ``verification_layers.rst``.

Pure APIs: ``TurbVisFluxPolicy.h``.
