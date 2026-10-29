Limited face reconstruction
===========================

After cell limiters are known, left/right face states are built by a limited
linear reconstruction. This sits on the residual path between gradient /
limiter evaluation and the numerical flux.

Pipeline
--------

.. code-block:: text

   qf = cell value (GetQlQr)
        |
        +--> phi = BlendFaceLimiter(phi_cell, phi_neighbor, both-sides?)
        +--> qf  = ApplyLimitedFaceReconstruction(qf, phi, g, x_f - x_c)
        |
        +-- if positivity check fails (weighted path):
              qf <- WeightedFaceAverage with InverseDistanceFaceWeights

Formula
-------

.. math::

   q_f = q_c + \phi\,\nabla q_c\cdot(\mathbf{x}_f - \mathbf{x}_c)

``phi`` is either the owner-cell limiter or
:math:`\min(\phi_L,\phi_R)` when face-two-sided blending is enabled
(``ireconface != 0``).

Viscous gradient correction (shared policy)
------------------------------------------

.. math::

   \nabla q \leftarrow \nabla q +
   \Big(\frac{q_R-q_L}{|LR|} - \nabla q\cdot\mathbf{u}\Big)\mathbf{u}

Pure APIs: ``ApplyLimitedFaceReconstruction``, ``BlendFaceLimiter``
(``ULimiterPolicy.h``); ``CorrectFaceGradient`` (``UGradPolicy.h``).
Tests: ``limiter_policy_test``, ``grad_policy_test``.
