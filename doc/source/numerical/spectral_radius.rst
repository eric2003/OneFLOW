Face spectral radius and local time step
========================================

Local time stepping uses face inviscid and viscous spectral radii accumulated
to cells. The pure formulas below are independent of mesh traversal.

Inviscid face spectral radius
-----------------------------

.. math::

   \lambda_{\mathrm{inv}}
   = \tfrac12 A \big(|v_n| + c\big)

with face-midpoint relative normal velocity and ideal-gas sound speed.

Viscous estimates
-----------------

Model 1 (momentum/energy diffusivity max):

.. math::

   \lambda_{\mathrm{vis}}
   = A^2 \cdot
   \frac{2\,\max\big(\tfrac43(\mu_l+\mu_t),\,
   \gamma(\mu_l/Pr_l+\mu_t/Pr_t)\big)}
   {Re\,\rho}

Model 2 uses face-to-cell distance; model 3 is a simplified
:math:`A^2/(Re\rho)` form.

Cell time step
--------------

.. math::

   \Delta t_{\mathrm{inv}} = \mathrm{CFL}\,\frac{V}{\lambda_{\mathrm{inv}}},
   \qquad
   \Delta t \leftarrow
   \frac{\Delta t_{\mathrm{inv}}\,\Delta t_{\mathrm{vis}}}
   {\Delta t_{\mathrm{inv}}+\Delta t_{\mathrm{vis}}}

References
----------

* Blazek, J. *Computational Fluid Dynamics: Principles and Applications*
  (spectral-radius estimates for viscous CFL).
* Toro, E. F. *Riemann Solvers and Numerical Methods for Fluid Dynamics*
  (wave speeds entering inviscid spectral radius).
* OpenFOAM Programmer's / User Guide (local time-stepping practice; compare
  implementation choices, not copy code).

Pure APIs: ``TimeStepPolicy.h``. Tests: ``timestep_policy_test``.
