Viscous conductive heat flux
============================

Face energy residual includes a conductive contribution aligned with the
temperature normal gradient.

Effective conductivity
----------------------

.. math::

   k = \Big(\mu_l / Pr_l + \mu_t / Pr_t\Big) c_p

Normal heat flux
----------------

.. math::

   q_n = \mathbf{q}\cdot\mathbf{n} + k\,\partial_n T

With Cartesian heat-vector components currently zero in the laminar path,
the dominant term is :math:`k\partial_n T`.

Pure APIs: ``EffectiveThermalConductivity``, ``ConductiveHeatFluxNormal``,
``FaceNormalDot`` (``NsVisFluxPolicy.h``). Tests: ``vis_flux_policy_test``.
