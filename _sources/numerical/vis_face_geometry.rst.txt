Viscous face geometry contracts
===============================

Shares inverse-distance face weights with Green-Gauss (``UGradPolicy``).
Viscous reconstruction (``VisGradGeom``) no longer re-implements the weight
formula; it calls the same pure function.

Pipeline (geometry side)
------------------------

.. code-block:: text

   cell/face centers
        |
        +--> InverseDistanceFaceWeights -> fw1, fw2  (blend cell grads)
        +--> FaceSkewAngleDegrees       -> angle1/2 (skew validity gates)
        +--> UnitizeDisplacement        -> unit L->R and 1/|LR|

Skew angle
----------

.. math::

   \theta^{\circ} = \frac{180}{\pi}
   \arcsin\Big(\mathrm{clamp}\big(\tfrac{p}{|\mathbf{d}|}, -1, 1\big)\Big)

where :math:`p` is the signed projection length used by the legacy vis path
(``d1``/``d2`` along the face normal).

Unit cell-to-cell direction
---------------------------

.. math::

   \mathbf{u}_{LR} = \frac{\mathbf{x}_{cR}-\mathbf{x}_{cL}}{|\mathbf{x}_{cR}-\mathbf{x}_{cL}|},
   \quad
   o_{ds} = \frac{1}{|\mathbf{x}_{cR}-\mathbf{x}_{cL}|}

Pure APIs live in ``UGradPolicy.h`` (``@doc-algorithm vis-face-geometry`` and
``green-gauss-gradient``). Tests: ``grad_policy_test``.
