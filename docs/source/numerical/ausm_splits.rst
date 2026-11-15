AUSM Mach and pressure splits
=============================

AUSM-family schemes split the face Mach number and pressure with low-order
polynomials in the subsonic range and upwind steps outside the sonic point.
These splits are pure and independent of mesh state.

Mach splits
-----------

With sense :math:`s=\pm 1` (left/right):

.. math::

   M^{\pm}_{(1)} = \tfrac{1}{2}(M \pm |M|)

.. math::

   M^{\pm}_{(2)} = \pm \tfrac{1}{4}(M \pm 1)^2

.. math::

   M^{\pm}_{(4)} =
   \begin{cases}
   M^{\pm}_{(1)} & |M| \ge 1 \\
   M^{\pm}_{(2)}\big(1 \mp 16\beta M^{\mp}_{(2)}\big) & |M| < 1
   \end{cases}

Pressure split
--------------

.. math::

   P^{\pm}_{(5)} =
   \begin{cases}
   \tfrac{1}{2}\big(1 \pm \mathrm{sign}(M)\big) & |M| \ge 1 \\
   M^{\pm}_{(2)}\big((2\mp M) \mp 16\alpha M M^{\mp}_{(2)}\big) & |M| < 1
   \end{cases}

Pure APIs: ``FMSplit1``, ``FMSplit2``, ``FMSplit4``, ``FPSplit5``
(``NsInvFluxPolicy.h``). Tests: ``inv_flux_policy_test``.
