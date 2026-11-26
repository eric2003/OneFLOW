Viscous face spectral radius (model 1)
======================================

Purpose
-------

Model 1 estimates the viscous contribution to the local spectral radius for
explicit or semi-implicit time-step control. It is a Blazek-style engineering
bound, not an exact Jacobian eigenvalue.

Formula
-------

With laminar viscosity :math:`\mu`, turbulent viscosity :math:`\mu_t`, face
area :math:`A`, density :math:`\rho`, and Reynolds number :math:`\mathrm{Re}`:

.. math::

   c_1=\frac{4}{3}(\mu+\mu_t),\qquad
   c_2=\gamma\left(\frac{\mu}{\mathrm{Pr}_l}+\frac{\mu_t}{\mathrm{Pr}_t}\right)

   \lambda_{\mathrm{vis}}=A^{2}\cdot\frac{2\max(c_1,c_2)}{\mathrm{Re}\,\rho}

OneFLOW stores inverse Prandtl numbers, so
:math:`c_2=\gamma(\mu\cdot\mathrm{invPr}_l+\mu_t\cdot\mathrm{invPr}_t)`.

Derivation
----------

Diffusion model
^^^^^^^^^^^^^^^

A one-dimensional diffusion equation :math:`u_t=\nu u_{xx}` damps a Fourier
mode with rate :math:`\nu k^{2}`. On a finite-volume cell, the mesh wave
number is estimated from the face-to-volume scale:

.. math::

   k \sim \frac{A}{V}
   \quad\Rightarrow\quad
   \text{decay rate} \sim \nu\left(\frac{A}{V}\right)^{2}.

For time-step control the viscous spectral contribution is therefore taken
proportional to :math:`A^{2}\nu_{\mathrm{eff}}`, with density and Reynolds
scaling so that the result has the same dimensions as the inviscid face
spectral radius and can enter

.. math::

   \Delta t \sim \frac{\mathrm{CFL}\cdot V}{\sum(\lambda_{\mathrm{inv}}+\lambda_{\mathrm{vis}})}.

Momentum channel
^^^^^^^^^^^^^^^^

The effective dynamic viscosity is :math:`\mu+\mu_t`. A three-dimensional
engineering factor :math:`4/3` accounts for longitudinal / bulk-viscous
effects in the spectral estimate:

.. math::

   c_1=\frac{4}{3}(\mu+\mu_t).

Energy channel
^^^^^^^^^^^^^^

Heat conduction uses laminar and turbulent Prandtl numbers. With an ideal-gas
factor :math:`\gamma` the energy diffusivity scale is

.. math::

   c_2=\gamma\left(\frac{\mu}{\mathrm{Pr}_l}+\frac{\mu_t}{\mathrm{Pr}_t}\right).

Assembly
^^^^^^^^

The tighter (larger) of the two channels is retained, scaled by
:math:`2/(\mathrm{Re}\,\rho)` and by :math:`A^{2}`:

.. math::

   \lambda_{\mathrm{vis}}
   =A^{2}\cdot\frac{2\max(c_1,c_2)}{\mathrm{Re}\,\rho}.

The factor :math:`2` is an engineering constant of the estimate, not a unique
eigenvalue coefficient from a single continuous operator.

Relation to inviscid spectrum
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The inviscid face contribution remains

.. math::

   \lambda_{\mathrm{inv}}=\frac{1}{2}A\bigl(|v_n|+c\bigr)

(``InviscidFaceSpectralRadius``). Cell time steps may combine inviscid and
viscous limits through ``CellTimeStepFromSpectrum`` and
``CombineInvAndVisTimeStep``.

Source map
----------

* ``ViscousFaceSpectralRadiusModel1`` in ``TimeStepPolicy``
* Related: ``InviscidFaceSpectralRadius``, ``CellTimeStepFromSpectrum``
* Tests: ``tests/solver/timestep_policy_test.cpp``
  (``ViscousFaceSpectralRadiusModel1CurveMatchesClosedForm``)
* Theory table: ``tools/verification/data/viscous_spectral_radius_model1_theory.csv``
* Contract Workbench id: ``viscous_spectral_radius_model1``

.. image:: figures/viscous_spectral_radius_model1_theory.png
   :alt: Viscous spectral radius model1 vs mut
   :width: 520

What this proves / does not prove
---------------------------------

**Proves:** model-1 algebra matches the pure API on a sample :math:`\mu_t` grid;
scaling with :math:`A^{2}` holds in the unit test.

**Does not prove:** exact Jacobian spectral radius; optimality of CFL on a
case mesh; models 2 and 3.

References
----------

* J. Blazek, *Computational Fluid Dynamics: Principles and Applications*,
  Elsevier (2nd/3rd ed.), temporal discretisation / maximum time step:
  viscous spectral radius with
  \(\maxigl(	frac{4}{3}(\mu_L+\mu_T),\gamma(\mu_L/\mathrm{Pr}_L+\mu_T/\mathrm{Pr}_T)igr)\)
  and surface-element scaling. Primary structural source for model 1.
* Overview and comparative notes: :doc:`viscous_spectral_radius`.
* SU2 / AIAA papers: edge-based \(\lambda_{\mathrm{vis}}\propto\mu S^{2}/ho\)
  for comparison of scaling only.
* OneFLOW intermediate curve:
  ``tools/verification/data/viscous_spectral_radius_model1_theory.csv``.
