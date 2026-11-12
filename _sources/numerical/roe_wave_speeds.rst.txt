Roe wave speeds and entropy fix
===============================

After the Roe average state is known, inviscid schemes need the face-normal
wave speeds and (optionally) a Harten entropy floor.

Formulas
--------

Total enthalpy:

.. math::

   H = h_{\mathrm{int}} + \tfrac{1}{2}|\mathbf{v}|^2

Ideal-gas sound speed from total enthalpy:

.. math::

   c = \sqrt{\left|( \gamma-1)\big(H - \tfrac{1}{2}|\mathbf{v}|^2\big)\right|}

Absolute Roe eigenvalues:

.. math::

   a_1 = |v_n|,\quad a_2 = |v_n+c|,\quad a_3 = |v_n-c|

Harten floor (when enabled): let :math:`a_{\max}=\max(a_2,a_3)`.

.. math::

   \tilde{a}_1 = \max(c_1 a_{\max}, a_1),\quad
   \tilde{a}_{2,3} = \max(c_2 a_{\max}, a_{2,3})

Pipeline
--------

.. code-block:: text

   RoeAverageState
        |
        +--> TotalEnthalpy (left/right)
        +--> IdealGasSoundSpeed / c2 at mid state
        +--> vn_rel at mid; RoeAbsoluteEigenvalues
        +--> ApplyHartenEntropyFix (ctrl.ieigenfix)
        +--> dissipation terms using meig*

Pure APIs in ``NsInvFluxPolicy.h``. Tests: ``inv_flux_policy_test``.
