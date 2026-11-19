Residual face-to-cell assembly
==============================

After inviscid/viscous face fluxes are computed, residuals accumulate with a
fixed orientation. This is the last pure contract on the spatial residual path
before linear algebra / time integration.

Convention
----------

.. code-block:: text

   face flux F oriented left -> right
        |
        +-- boundary face:  R_left  -= F
        |                   R_right unchanged
        |
        +-- interior face:  R_left  -= F
                            R_right += F

Formula
-------

.. math::

   R_L \leftarrow R_L - F_f,\qquad
   R_R \leftarrow R_R + F_f \quad (\text{interior only})

The same rule is used for inviscid and viscous residual fields
(``AddF2CField``).

Pure API: ``DepositOrientedFaceFlux`` in ``UResidualPolicy.h``
(``@doc-algorithm residual-assembly``). Tests: ``residual_deposit_test``.
