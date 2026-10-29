Euler convective face flux
==========================

Physical (non-dissipative) Euler flux on a face appears in Roe, AUSM-family and
related schemes before upwind dissipation is added. The formula is shared; only
the dissipation terms differ by scheme.

Pipeline role
-------------

.. code-block:: text

   left/right primitives
        |
        +--> vn_rel = n·v - v_fn
        +--> F_phys = EulerConvectiveFlux(rho, v, p, h, vn_rel, n, v_fn)
        |
        +--> scheme-specific dissipation (Roe matrix, AUSM splits, ...)
        +--> DepositOrientedFaceFlux into residual

Formula
-------

.. math::

   v_n^{\mathrm{rel}} = \mathbf{n}\cdot\mathbf{v} - v_{fn}

.. math::

   \mathbf{F} = \dot{m}
   \begin{bmatrix} 1 \\ u \\ v \\ w \\ h \end{bmatrix}
   + p
   \begin{bmatrix} 0 \\ n_x \\ n_y \\ n_z \\ v_{fn} \end{bmatrix},
   \quad \dot{m} = \rho\, v_n^{\mathrm{rel}}

Pure APIs: ``RelativeFaceNormalVelocity``, ``EulerConvectiveFlux``
(``NsInvFluxPolicy.h``). Tests: ``inv_flux_policy_test``.
