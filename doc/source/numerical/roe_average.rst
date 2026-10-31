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

LUSGS face-state closure
------------------------

The face state prepared by ``UNsLusgs::PrepareDataFacePrim()`` has a
solver-specific use. In ``NsLusgs::GetStandardFluxIncrement()``, the
prepared face state is consumed by ``CalcFaceEigenValue()`` to estimate the
face spectral radius. It is not used there to recompute the steady physical
flux residual. The resulting eigenvalue participates in the implicit LUSGS
correction. The current ``PrepareData()`` path subsequently assigns
``nslu.gama`` from the right cell, and ``CalcFaceEigenValue()`` uses that
value for the face sound speed. Thus the pressure closure uses left-cell
``gamma`` while this eigenvalue calculation uses right-cell ``gamma``. This
is existing behavior worth preserving in compatibility tests and auditing
separately for variable-gamma flows; it should not be silently normalized as
part of an unrelated refactor.

The current compressible unstructured LUSGS implementation first computes the
Roe-weighted total enthalpy and velocity, then closes pressure using the
left-cell heat-capacity ratio:

.. math::

   p_m^{H} = \rho_m
       \left(H_m - \frac{1}{2}|\mathbf{u}_m|^2\right)
       \frac{\gamma_L-1}{\gamma_L}

A weighted static-pressure alternative is:

.. math::

   p_m^{p} = \frac{p_L + r p_R}{1+r}

These expressions are not algebraically equivalent in general. For example,
with left state ``(rho=1, p=1, u=0)``, right state
``(rho=4, p=4, u=3)``, and ``gamma=1.4``, the current enthalpy closure
gives `18/7`, while the weighted static pressure gives `3`.

The difference changes the face sound speed and therefore the implicit
spectral-radius estimate. Since this face state is used in the LUSGS
linearization/preconditioning path, a different intermediate residual history
does not by itself prove that the underlying steady discretization has
changed or that one closure is mathematically wrong. If both iterations
converge to the same root of the same residual equations, their final
solutions may agree even when their iteration paths differ.

This is a hypothesis to verify, not an assertion that both choices are
equally robust. A meaningful A/B study should keep the mesh, physical
residual, boundary conditions, initialization, and convergence criteria
fixed, then compare:

* residual histories and convergence rates, not just one intermediate file;
* whether both runs reach the same converged residual level;
* conservative variables and derived quantities such as lift and drag;
* conservation balance and sensitivity to CFL, initialization, and mesh.

Regression result files serve a different purpose: they pin the established
numerical behavior for compatibility. The unit test for
``RoePressureFromTotalEnthalpy`` verifies the selected algebraic closure;
it does not establish that a weighted-static-pressure closure is universally
invalid. Any intentional change of closure should therefore be evaluated as
a numerical-method change, with convergence and solution comparisons, rather
than silently redefining the compatibility baseline.
