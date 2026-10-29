Roe average state
=================

Roe-type inviscid fluxes evaluate dissipation at an intermediate state.
The density is a geometric mean; other scalars use the Roe density-weighted
average.

Formulas
--------

.. math::

   r = \sqrt{\frac{\rho_R}{\rho_L}},\qquad
   c = \frac{1}{1+r}

.. math::

   \rho_m = \sqrt{\rho_L\rho_R},\qquad
   \phi_m = (\phi_L + \phi_R\, r)\, c
   \quad (\phi \in \{u,v,w,p,h,\gamma\})

Pipeline role
-------------

.. code-block:: text

   left/right primitives + total enthalpies
        |
        +--> RoeAverageState -> (rho,u,v,w,p,h,gamma)_m
        +--> sound speed / eigenvalues at mid state
        +--> physical flux + Roe dissipation
        +--> residual deposit

Pure APIs: ``RoeAverageWeights``, ``RoeAverageComponent``, ``RoeAverageState``
(``NsInvFluxPolicy.h``). Tests: ``inv_flux_policy_test``.
