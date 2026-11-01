SST eddy viscosity and F2
=========================

.. math::

   \mathrm{arg}_2
   = \max\Big(
   \frac{2\sqrt{k}}{\beta^* d\omega},\,
   \frac{500\nu}{\omega d^2 Re}
   \Big),
   \qquad
   F_2 = \tanh(\mathrm{arg}_2^2)

.. math::

   \mu_t
   = \frac{\rho k}
   {\max(\omega,\,R F_2/a_1)}\,Re

Here :math:`R` is the rate magnitude selected by the SST production option:
the strain-rate magnitude when ``iprod_sst == 1``, otherwise the
vorticity magnitude. For velocity gradient :math:`\nabla\mathbf{u}`, the
strain invariant is :math:`|S|=\sqrt{2 S_{ij}S_{ij}}`, with
:math:`S_{ij}=\tfrac12(\partial_j u_i+\partial_i u_j)`. The alternative
vorticity magnitude is :math:`|\Omega|=|\nabla\times\mathbf{u}|`. Both are
now computed by separately tested pure policy functions, and the selected
quantity is passed unchanged to the Bradshaw limiter. The solver applies an
optional laminar-viscosity
ratio cap after this expression, implemented by the shared
``TurbCappedEddyViscosity()`` policy used by both SA and SST.

Intermediate checks cover both denominator branches: the :math:`\omega`
limit when :math:`R F_2/a_1 \leq \omega`, and the strain/vorticity
limit when :math:`R F_2/a_1 > \omega`. The curve test also checks
:math:`F_2(0)=0` and :math:`F_2\to 1` for large arg.

References
----------

* Menter SST; NASA TMR SST.
* See ``sa_eddy_viscosity.rst`` and ``verification_layers.rst``.

Pure APIs: ``SstStrainRateMagnitude``, ``SstVorticityMagnitude``,
``SstF2Arg``, ``SstF2FromArg``, ``SstEddyViscosity``.

Intermediate verification
-------------------------

Closed-form checks of :math:`F_2(\mathrm{arg}_2)=\tanh(\mathrm{arg}_2^2)`
live in ``turb_sst_eddy_viscosity_test`` and theory table
``tools/verification/data/sst_f2_arg2_theory.csv`` (Contract Workbench Plot).

.. image:: figures/sst_f2_arg2_theory.png
   :alt: SST F2(arg2) theory curve
   :width: 480
