LUSGS approximate flux Jacobian increment
=========================================

The lower/upper sweep builds a face flux increment
:math:`\Delta F \approx (\partial F/\partial U)\Delta U` with a
spectral split of the Roe eigenvalues.

Split eigenvalues
-----------------

.. math::

   \lambda_i^{\pm}
   = \tfrac12\big(\lambda_i \pm \mathrm{sign}\,|\lambda|_{\mathrm{face},i}\big)

Linearization weights match the Roe dissipation structure
(:math:`x_1,x_2` analogous to :math:`\xi_1,\xi_2`).

Five-equation write-back is area-scaled and shares
``RoeDissipationCouplingTerms`` / ``RoeTotalEnthalpyJump``.

References
----------

* Yoon, S. and Jameson, A. Lower-upper symmetric Gauss-Seidel methods.
* Blazek, J. *Computational Fluid Dynamics: Principles and Applications*
  (implicit residual linearization).
* See also ``references.rst`` and ``roe_dissipation.rst``.

Pure APIs: ``NsLusgsPolicy.h``. Tests: ``lusgs_policy_test``.
