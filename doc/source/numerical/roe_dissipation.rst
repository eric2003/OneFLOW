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

Total-enthalpy jump
-------------------

.. math::

   \Delta H = -(\gamma-1)\big(\mathbf{v}\cdot\Delta(\rho\mathbf{v}) - \Delta(\rho e)\big)
   + \tfrac12(\gamma-1)|\mathbf{v}|^2\,\Delta\rho

Coupling terms and five-equation write-back
------------------------------------------

.. math::

   \mathrm{term}_1 = \Delta H\,\xi_1 + d_c\,\xi_2,\qquad
   \mathrm{term}_2 = c^2 d_c\,\xi_1 + \Delta H\,\xi_2

.. math::

   \mathbf{F} \leftarrow \mathbf{F}
   - \big(\tilde{a}_1\Delta\mathbf{q} - \mathbf{r}(\mathrm{term}_1,\mathrm{term}_2)\big)

with the standard Roe right-eigenvector combination encoded in
``ApplyRoeDissipationFiveEq``.
