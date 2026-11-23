Sutherland viscosity law
========================

Laminar dynamic viscosity is evaluated with the Sutherland formula, both
in dimensional form (reference tables) and nondimensional form (solver
working variables).

Nondimensional form
-------------------

With :math:`T` and :math:`t_s` already scaled by the reference temperature:

.. math::

   \frac{\mu}{\mu_{\mathrm{ref}}}
   = T^{3/2}\frac{1+t_s}{T+t_s}

Dimensional form
----------------

.. math::

   \mu(T)
   = \mu_0
   \left(\frac{T_0+T_s}{T+T_s}\right)
   \left(\frac{T}{T_0}\right)^{3/2}

Default air constants in ``Sutherland::Init``:
:math:`T_0=273.16\,\mathrm{K}`, :math:`T_s=110.4\,\mathrm{K}`,
:math:`\mu_0=1.715\times 10^{-5}\,\mathrm{Pa\,s}`.

Pure APIs: ``SutherlandViscosityNondim``, ``SutherlandViscosityDim``
(``Sutherland.h``). Class methods delegate to these. Tests:
``sutherland_policy_test``.
