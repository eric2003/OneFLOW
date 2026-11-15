Limiter face contracts
======================

ROI priority 1 (solver / numerical): the slope limiter sits on every
reconstructed face value. Pure scalar policies and neighbor-range construction
are already isolated; this page locks the **face geometry** steps shared by
Barth and Venkat.

Pipeline
--------

.. code-block:: text

   [neighbor offset bounds] --> [per face]
        |                          |
        |                          +--> dq_f = g . (x_f - x_c)
        |                          +--> (Barth only) s = orient * (n . d) / (|d|+eps)
        |                          +--> phi_face = Barth/Venkat(dq_f, ...)
        |                          +--> phi_cell = min(phi_cell, phi_face)
        v
   global eps for Venkat from offset range (once per equation)

Formulas
--------

Unlimited face increment:

.. math::

   \Delta q_f = \nabla q_c \cdot (\mathbf{x}_f - \mathbf{x}_c)

Barth side direction factor:

.. math::

   s = \sigma \frac{\mathbf{n}\cdot\mathbf{d}}{|\mathbf{d}| + \varepsilon},
   \quad \sigma = +1 \text{ (left)},\ -1 \text{ (right)}

Venkat regularization (global, once per field):

.. math::

   \varepsilon_{\mathrm{reg}} = \big(c \cdot \max(\Delta_{\max}-\Delta_{\min}, 1)\big)^2 + \delta

where :math:`\Delta` are global neighbor **offset** ranges.

Pure entry points (``ULimiterPolicy.h``)
----------------------------------------

* ``ReconstructFaceIncrement``
* ``FaceSideDirectionDot``
* ``TightenLimiter``
* ``BarthFunction`` / ``VencatFunction`` / ``VencatRegularizationEpsilon``
* ``BuildCellNeighborOffsetBounds`` (``ULimiterRange.h``)

Header tags ``@doc-algorithm slope-limiter-face`` mark the face stages for
optional extraction. Tests: ``limiter_policy_test``, ``limiter_range_test``.
