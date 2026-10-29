Green-Gauss gradient face contracts
===================================

Upstream of the slope limiter on the residual path. Cell gradients are built by
accumulating face contributions; the **face value** and **area-weighted normal**
are pure and testable without mesh globals.

Pipeline
--------

.. code-block:: text

   for each face:
     q_f  <- arithmetic average   OR   inverse-distance weights
     S    <- n * A
     left cell  += S * q_f
     interior right cell -= S * q_f
   for each physical cell:
     grad *= 1 / volume
   boundary ghosts copy owner gradient

Inverse-distance weights
------------------------

.. math::

   d_L = |\mathbf{x}_f - \mathbf{x}_{cL}|,\quad
   d_R = |\mathbf{x}_f - \mathbf{x}_{cR}|

.. math::

   c_L = \frac{d_R}{d_L + d_R + \varepsilon},\quad
   c_R = \frac{d_L}{d_L + d_R + \varepsilon},\quad
   q_f = c_L q_L + c_R q_R

Closer cell gets larger weight. Equal distances recover the arithmetic average.

Pure entry points (``UGradPolicy.h``)
------------------------------------

* ``ArithmeticFaceAverage``
* ``InverseDistanceFaceWeights``
* ``WeightedFaceAverage``
* ``AreaWeightedFaceNormal``

Tagged ``@doc-algorithm green-gauss-gradient``. Tests: ``grad_policy_test``.
