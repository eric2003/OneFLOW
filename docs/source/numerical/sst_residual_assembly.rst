SST residual assembly (Menter linearization)
============================================

After production and dissipation are known, the two-equation residual and
implicit diagonal terms are assembled as pure combinations.

.. math::

   S_k = P_k - \varepsilon,
   \qquad
   S_\omega = P_\omega - D_\omega + CD_{k\omega}

.. math::

   d_k = -2 f_{\beta^*}\beta^*\omega,
   \qquad
   d_\omega = -2 f_\beta\beta\omega - |CD|/(\rho\omega)

Negative parts :math:`\tfrac12(S-|S|)` feed the implicit sweep.

References
----------

* Menter, F. R. SST model residual linearization practice.
* NASA TMR — SST.
* See ``verification_layers.rst``.

Pure APIs: ``SstSourceK``, ``SstSourceW``, ``SstDiagonalK``,
``SstDiagonalW``, ``SstNegativeSourcePart``.
