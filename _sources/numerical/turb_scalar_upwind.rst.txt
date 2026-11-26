Turbulence scalar convective upwind flux
========================================

Turbulence working variables (SA :math:`\tilde\nu`, SST :math:`k,\omega`, ...)
use a face-normal upwind convective flux with optional density weighting.

Relative normal velocity
------------------------

.. math::

   v_n = \mathbf{n}\cdot\mathbf{v} - v_{fn}

Split
-----

.. math::

   v_n^+ = \tfrac12(v_n+|v_n|),\qquad
   v_n^- = \tfrac12(v_n-|v_n|)

Density weight
--------------

.. math::

   w = \alpha + (1-\alpha)\rho

with :math:`\alpha=1` for SA (flux of :math:`\tilde\nu`) and
:math:`\alpha=0` for density-weighted two-equation models.

Scalar flux
-----------

.. math::

   F = (v_n^+ w_L)\phi_L + (v_n^- w_R)\phi_R

References
----------

* Spalart, P. R. and Allmaras, S. R. A one-equation turbulence model.
* Menter, F. R. Two-equation eddy-viscosity turbulence models (SST).
* See ``references.rst``.

Pure APIs: ``TurbInvFluxPolicy.h``. Tests: ``turb_inv_flux_policy_test``.
