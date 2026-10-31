Turbulence LUSGS scalar sweep kernels
=====================================

Turbulence implicit sweeps treat each equation as a scalar system with
face matrix contributions and a cell spectral radius.

.. math::

   \Delta F = M\,\Delta q,
   \qquad
   \Delta q = r / \lambda

Sweep convergence accumulates :math:`(\Delta q-\Delta q_0)^2`.

References
----------

* Yoon–Jameson LUSGS; turbulence scalar implicit practice.
* See ``lusgs_flux_increment.rst`` (NS) and ``turb_spectrum.rst``.

Pure APIs: ``TurbLusgsPolicy.h``. Tests: ``turb_lusgs_policy_test``.
