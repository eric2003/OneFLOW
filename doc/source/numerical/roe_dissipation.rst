Roe dissipation coefficients
============================

After physical fluxes and absolute eigenvalues (with optional Harten fix)
are known, the Roe scheme subtracts a dissipative correction built from
wave coefficients and jumps in conserved variables.

Wave coefficients
-----------------

.. math::

   \xi_1 = \frac{2\tilde{a}_1 - \tilde{a}_2 - \tilde{a}_3}{2 c^2},\qquad
   \xi_2 = \frac{\tilde{a}_2 - \tilde{a}_3}{2 c}

where :math:`\tilde{a}_i` are the entropy-fixed absolute eigenvalues.

Continuity / normal mass jump used in the energy coupling:

.. math::

   d_c = v_n\,\Delta\rho - \mathbf{n}\cdot\Delta(\rho\mathbf{v})

Conservative to primitive (5-eq)
--------------------------------

.. math::

   \rho = |q_\rho|,\quad
   \mathbf{v} = \mathbf{q}_{\rho v}/\rho,\quad
   p = (\gamma-1)\big(q_E - \tfrac12\rho|\mathbf{v}|^2\big)

Pure APIs: ``RoeDissipationWaveCoefficients``, ``RoeContinuityNormalJump``,
``IdealGasConservativeToPrim`` (``NsInvFluxPolicy.h``).
Tests: ``inv_flux_policy_test``.
