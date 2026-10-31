SST eddy viscosity and F2
=========================

.. math::

   \mathrm{arg}_2
   = \max\Big(
   \frac{2\sqrt{k}}{\beta^* d\omega},\,
   \frac{500\nu}{\omega d^2 Re}
   \Big),
   \qquad
   F_2 = \tanh(\mathrm{arg}_2^2)

.. math::

   \mu_t
   = \frac{\rho k}
   {\max(\omega,\,\Omega F_2/a_1)}\,Re

with an optional laminar-ratio cap.

Intermediate checks: \(F_2(0)=0\), \(F_2\to 1\) for large arg.

References
----------

* Menter SST; NASA TMR SST.
* See ``sa_eddy_viscosity.rst`` and ``verification_layers.rst``.

Pure APIs: ``SstF2Arg``, ``SstF2FromArg``, ``SstEddyViscosity``.
