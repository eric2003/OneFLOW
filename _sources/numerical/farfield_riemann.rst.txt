Subsonic far-field Riemann invariants
=====================================

For subsonic far-field faces, one-dimensional acoustic Riemann invariants
along the face normal close the boundary state.

Invariants
----------

.. math::

   R^+ = v_n^{\mathrm{in}} + \frac{2 c_{\mathrm{in}}}{\gamma-1},
   \qquad
   R^- = v_n^{\mathrm{ref}} - \frac{2 c_{\mathrm{ref}}}{\gamma-1}

.. math::

   v_{n,b} = \tfrac12(R^+ + R^-),
   \qquad
   c_b = \tfrac14(R^+ - R^-)(\gamma-1)

Tangential velocity is taken from the upwind side (exit: interior; inlet:
freestream). Density follows isentropic relation
:math:`\rho = (c_b^2/(s\gamma))^{1/(\gamma-1)}` and
:math:`p = c_b^2\rho/\gamma`.

References
----------

* Blazek, J. *Computational Fluid Dynamics: Principles and Applications*.
* Toro, E. F. *Riemann Solvers and Numerical Methods for Fluid Dynamics*.
* See ``boundary_velocity.rst`` and ``references.rst``.

Pure APIs: ``NsBcPolicy.h``. Tests: ``bc_policy_test``.
