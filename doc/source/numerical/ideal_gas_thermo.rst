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
