Spalart-Allmaras closure functions
==================================

The SA one-equation model evaluates a chain of nondimensional closure
functions from the working variable :math:`\tilde\nu` and the local
strain/vorticity magnitude.

Core definitions
----------------

.. math::

   \chi = \frac{\tilde\nu}{\nu},
   \qquad
   f_{v1} = \frac{\chi^3}{\chi^3+c_{v1}^3},
   \qquad
   f_{v2} = 1 - \frac{\chi}{1+\chi f_{v1}}

Modified strain
---------------

With :math:`r_s = \tilde\nu/(\kappa^2 d^2)` (including Reynolds scaling
in the solver) and :math:`\bar S = r_s f_{v2}`:

.. math::

   \tilde S =
   \begin{cases}
   \Omega + \bar S & \bar S \ge -c_2\Omega \\
   \Omega\big(1 + \frac{c_2^2\Omega + c_3\bar S}
   {(c_3-2c_2)\Omega - \bar S}\big) & \text{otherwise}
   \end{cases}

Wall destruction function
-------------------------

.. math::

   r = \min(r_s/\tilde S,\,10),
   \quad
   g = r + c_{w2}(r^6-r),
   \quad
   f_w = g\Big(\frac{1+c_{w3}^6}{g^6+c_{w3}^6}\Big)^{1/6}

Optional trip function :math:`f_{t2}=c_{t3}\exp(-c_{t4}\chi^2)`.

References
----------

* Spalart, P. R. and Allmaras, S. R. (1992/1994). A one-equation
  turbulence model for aerodynamic flows.
* Allmaras, S. R., Johnson, F. T., Spalart, P. R. Modifications and
  clarifications for the implementation of the Spalart-Allmaras model.
* NASA Turbulence Modeling Resource (SA equations and constants).
* See ``references.rst`` for the citation policy and self-evolution loop.

Pure APIs: ``TurbSaPolicy.h``. Tests: ``turb_sa_policy_test``.
Wiring: ``TurbCom::CalcSrcSa``.
