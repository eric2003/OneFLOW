Isothermal wall density and temperature images
==============================================

Viscous isothermal walls set the face temperature to :math:`T_w` and build
ghost density/temperature by linear images subject to floors.

Equation of state
-----------------

.. math::

   \rho = \frac{p}{R\,T\,\overline{W}^{-1}}

Temperature image
-----------------

.. math::

   T' = \max\big(T_{\min},\, 2 T_w - T_{\mathrm{in}}\big)

Density image (first ghost layer)
---------------------------------

.. math::

   \rho' = \max\big(\tfrac12\rho_w,\, 2\rho_w - \rho_{\mathrm{in}}\big)

where :math:`\rho_w = p/(R T_w)`.

References
----------

* Blazek, J. *Computational Fluid Dynamics: Principles and Applications*
  (isothermal wall treatment).
* See ``boundary_velocity.rst`` and ``references.rst``.

Pure APIs: ``IdealGasDensityFromPressureTemperature``,
``MirrorScalarWithFloor``, ``MirrorDensityAboutWall`` (``NsBcPolicy.h``).
Tests: ``bc_policy_test``.
