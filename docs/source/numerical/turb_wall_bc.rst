Turbulence viscous-wall boundary values
=======================================

SA (one equation)
-----------------

Ghost working variable is the linear image of the interior so the
wall value is zero:

.. math::

   \tilde\nu_g = - \tilde\nu_i

SST / k-omega
-------------

Kinetic energy uses the same image; specific dissipation follows Menter:

.. math::

   \omega_w = \frac{60\mu_l}{\rho\beta_1 d^2 Re}

References
----------

* Menter SST wall boundary treatment.
* Spalart–Allmaras wall condition.
* NASA TMR.
* See ``verification_layers.rst``.

Pure APIs: ``TurbBcPolicy.h``. Tests: ``turb_bc_policy_test``.
