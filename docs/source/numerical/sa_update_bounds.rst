SA working-variable bounds after update
=======================================

After each explicit/implicit update, the SA working variable is limited.

Non-negativity
--------------

Negative values are replaced by a fallback (typically freestream
:math:`\tilde\nu`):

.. math::

   \tilde\nu \leftarrow
   \begin{cases}
   \tilde\nu_\infty & \tilde\nu < 0 \\
   \tilde\nu & \text{otherwise}
   \end{cases}

Eddy-viscosity ratio cap
------------------------

With :math:`f_{v1}(\chi)` evaluated at the current state:

.. math::

   \tilde\nu_{\max}
   = r_{\max}\,\frac{\mu_l}{\rho f_{v1}}

References
----------

* Spalart–Allmaras implementation practice; NASA TMR.
* See ``sa_eddy_viscosity.rst`` and ``verification_layers.rst``.

Pure APIs: ``SaWorkingVariableUpperBound``,
``SaNonNegativeWorkingVariable``.
