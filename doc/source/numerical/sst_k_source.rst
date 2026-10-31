SST k-equation production and dissipation
=========================================

Two-equation SST models share a pure production/dissipation kernel for
the turbulent kinetic energy equation.

Production variants
-------------------

* **Boussinesq** (full):
  :math:`P_k = \mu_t(S_{ij}:S_{ij}-\tfrac23(\nabla\cdot\mathbf{v})^2)/Re
  - \tfrac23\rho k\nabla\cdot\mathbf{v}`
* **SST-2003** strain-only: :math:`P_k = \mu_t S_{ij}:S_{ij}/Re`
* **SST-V** vorticity: replace :math:`S:S` by vorticity magnitude squared.

Dissipation
-----------

.. math::

   \varepsilon = f_{\beta^*}\beta^*\rho k\omega

or DES form :math:`f_{\beta^*}\rho k^{3/2}/L`. Production is limited by
:math:`P_k \leftarrow \min(P_k,\, p_{k,\mathrm{lim}}\varepsilon)`.

References
----------

* Menter, F. R. Two-equation eddy-viscosity turbulence models for
  engineering applications (SST).
* NASA Turbulence Modeling Resource — SST
  (https://turbmodels.larc.nasa.gov/sst.html).
* See ``sa_source.rst`` and ``references.rst`` for the citation /
  intermediate-verification policy.

Pure APIs: ``TurbSstPolicy.h``. Tests: ``turb_sst_policy_test``.
