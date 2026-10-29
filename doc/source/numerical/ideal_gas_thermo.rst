Ideal-gas thermodynamic contracts
=================================

Inviscid flux paths repeatedly convert between static enthalpy, total
enthalpy, internal energy and sound speed. For ``chemModel == 0`` these are
closed-form ideal-gas relations and are pure.

Formulas
--------

Static enthalpy:

.. math::

   h = \frac{\gamma}{\gamma-1}\frac{p}{\rho}

Specific total internal energy (used by ``PrimToQ``):

.. math::

   e = \frac{p}{\rho(\gamma-1)} + \tfrac{1}{2}|\mathbf{v}|^2

Sound speed from pressure:

.. math::

   c = \sqrt{\gamma p / \rho}

Relation to total enthalpy:

.. math::

   H = h + \tfrac{1}{2}|\mathbf{v}|^2
   = e + \frac{p}{\rho}

Pipeline role
-------------

.. code-block:: text

   primitives (rho, v, p)
        |
        +--> IdealGasStaticEnthalpy -> hint
        +--> TotalEnthalpy(hint, |v|^2) -> H
        +--> IdealGasSpecificInternalEnergy -> e  (PrimToQ)
        +--> IdealGasSoundSpeedFromPressure -> c  (Mach, AUSM, Van Leer)

Pure APIs in ``NsInvFluxPolicy.h``. Tests: ``inv_flux_policy_test``.

Conservative variables (5-eq)
-----------------------------

.. math::

   \mathbf{q} = \big(\rho,\;\rho u,\;\rho v,\;\rho w,\;\rho e\big),
   \quad
   e = \frac{p}{\rho(\gamma-1)} + \tfrac{1}{2}|\mathbf{v}|^2

Face Mach number:

.. math::

   M = v_n^{\mathrm{rel}} / c

``PrimToQ`` uses ``IdealGasPrimToConservative`` for the five core equations;
extra species equations remain ``rho * prim[i]`` in the caller.


Bilateral face sound speeds (AUSM-family)
-----------------------------------------

.. math::

   c_L = \sqrt{\gamma_L p_L/\rho_L},\quad
   c_R = \sqrt{\gamma_R p_R/\rho_R},\quad
   c_{\mathrm{avg}} = \tfrac12(c_L+c_R)

Some schemes replace :math:`c_{\mathrm{avg}}` with :math:`\max(c_L,c_R)`.
