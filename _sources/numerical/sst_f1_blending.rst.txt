SST blending function F1
========================

Menter SST switches between inner (Wilcox) and outer (k-epsilon-like)
constants with a continuous blending function :math:`F_1`.

Definitions
-----------

.. math::

   CD_{k\omega}
   = 2\rho\sigma_{\omega 2}
   (\nabla k\cdot\nabla\omega)/\omega

.. math::

   \mathrm{arg}_1
   = \min\Big(
   \max\Big(
   \frac{\sqrt{k}}{\beta^* d\omega},\,
   \frac{500\nu}{\omega d^2 Re}
   \Big),\,
   \frac{4\rho\sigma_{\omega 2} k}{CD_{k\omega} d^2}
   \Big)

.. math::

   F_1 = \tanh(\mathrm{arg}_1^4)

Intermediate verification
-------------------------

* Algebra: pure APIs for :math:`CD`, :math:`\mathrm{arg}_1`, :math:`F_1`.
* Curve: :math:`F_1(\mathrm{arg}_1)\to 0` as :math:`\mathrm{arg}_1\to 0`
  and :math:`F_1\to 1` for large :math:`\mathrm{arg}_1` (gtest).
* Later: wall-normal :math:`F_1(y)` profiles vs NASA TMR / reference codes
  (CSV/plot), still driven by these APIs.

References
----------

* Menter, F. R. Two-equation eddy-viscosity turbulence models (SST).
* NASA TMR — SST (https://turbmodels.larc.nasa.gov/sst.html).
* See ``sst_k_source.rst``, ``sst_omega_source.rst``, ``references.rst``.

Pure APIs: ``SstCrossDiffScalar``, ``SstCdkwWithFloor``,
``SstBlendingArg1``, ``SstF1FromArg1``, ``SstGradientDotProduct``.
