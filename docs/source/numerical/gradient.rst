Green-Gauss gradient contracts
==============================

Upstream of the slope limiter on the residual path. Face values, normals,
volume scaling, and boundary ghost copies are pure; only the face/cell
iteration still touches mesh globals (``ug``).

Pipeline
--------

.. code-block:: text

   zero cell gradients
        |
        v
   for each face:
        q_f  <- ArithmeticFaceAverage  OR  InverseDistance + WeightedFaceAverage
        S    <- AreaWeightedFaceNormal (n * A)
        left  AccumulateFaceGradient(+S * q_f)
        if interior:
          right AccumulateFaceGradient(-S * q_f)
        |
        v
   for each physical cell:
        ScaleGradientByVolume (grad /= V)
        |
        v
   for each boundary face:
        CopyGradient(owner -> ghost)

Formulas
--------

Arithmetic face value:

.. math::

   q_f = \frac{1}{2}(q_L + q_R)

Inverse-distance weights:

.. math::

   d_L = |\mathbf{x}_f - \mathbf{x}_{cL}|,\quad
   d_R = |\mathbf{x}_f - \mathbf{x}_{cR}|

.. math::

   c_L = \frac{d_R}{d_L + d_R + \varepsilon},\quad
   c_R = \frac{d_L}{d_L + d_R + \varepsilon},\quad
   q_f = c_L q_L + c_R q_R

Green-Gauss accumulation and volume scale:

.. math::

   \nabla q_c \leftarrow \nabla q_c + \sigma\,\mathbf{S}\,q_f,
   \quad \sigma = +1\ \mathrm{(left)},\ -1\ \mathrm{(right)}

.. math::

   \nabla q_c \leftarrow \frac{1}{V_c}\nabla q_c

Boundary ghosts copy the owner (left) gradient; they are not independent
unknowns in this scheme.

Pure entry points (``UGradPolicy.h``)
------------------------------------

* ``ArithmeticFaceAverage``
* ``InverseDistanceFaceWeights`` / ``WeightedFaceAverage``
* ``AreaWeightedFaceNormal``
* ``AccumulateFaceGradient``
* ``ScaleGradientByVolume``
* ``CopyGradient``

Tags: ``@doc-algorithm green-gauss-gradient`` with ``@doc-stage`` /
``@doc-formula`` on each pure API so an extractor can rebuild this page from
headers without a second algorithm write-up. Narrative, diagrams, and math
remain in ``doc/source`` (not only in ``//`` comments).

Tests: ``grad_policy_test``.
