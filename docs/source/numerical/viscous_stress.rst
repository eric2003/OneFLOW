Stokes viscous stress and face flux
===================================

After face velocity gradients are known, the laminar viscous flux is built
from the Stokes stress tensor and the heat-flux normal component.

Stress tensor
-------------

.. math::

   \tau_{ij}
   = \mu\Big(\partial_j u_i + \partial_i u_j\Big)
   - \tfrac{2}{3}\mu\,(\nabla\cdot\mathbf{v})\,\delta_{ij}

Momentum and energy face fluxes
-------------------------------

.. math::

   \mathbf{F}_{\mathrm{mom}} = \boldsymbol{\tau}\cdot\mathbf{n}

.. math::

   F_E = \mathbf{v}\cdot\mathbf{F}_{\mathrm{mom}} + q_n

Pipeline
--------

.. code-block:: text

   face velocity/temperature gradients
        |
        +--> StokesStressTensor (optional anisotropy add-on remains in caller)
        +--> ProjectStressMomentumFlux
        +--> ViscousEnergyFlux (+ heat normal)
        +--> DepositOrientedFaceFlux into residual

Pure APIs: ``StokesStressTensor``, ``ProjectStressMomentumFlux``,
``ViscousEnergyFlux`` (``NsVisFluxPolicy.h``). Tests: ``vis_flux_policy_test``.


Nondimensional face scaling
---------------------------

.. math::

   \mathrm{coeff} = -\frac{1}{Re}\,A

so that stored face viscous fluxes match the residual convention
``DepositOrientedFaceFlux``. Spectrum accumulation uses the same bilateral
deposit helper ``DepositFaceScalarToCells``.
