SST and transition update floors
================================

After the density-scaled increment

.. math::

   q \leftarrow q_0 + \Delta q / \rho

kinetic energy and specific dissipation are floored:

.. math::

   k \leftarrow \max(k, k_{\min}),
   \qquad
   \omega \leftarrow \max(\omega, \omega_{\min})

With gamma–Re_theta transition, intermittency is capped and Re_theta floored.

References
----------

* Menter SST implementation practice; transition model bounds.
* See ``sa_update_bounds.rst`` and ``verification_layers.rst``.

Pure APIs: ``SstDensityScaledUpdate``, ``SstFloorKineticEnergy``,
``SstFloorOmega``, ``SstCapIntermittency``, ``SstFloorReTheta``.
