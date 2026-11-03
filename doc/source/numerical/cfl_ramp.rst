CFL ramp and local time-step control
=====================================

Purpose
-------

The parameters ``cflst``, ``cfled``, and ``ncfl`` define a CFL
schedule for the unstructured Navier-Stokes time-step path. They do not
directly specify a physical time increment.

CFL schedule
------------

``Iteration::Init`` reads the three values from the runtime parameter
database. ``UTimeStep::CalcTimeStep`` calls ``CalcCfl()`` before selecting
the CPU or supported HIP path. ``TimeStep::CalcCfl`` uses
``outerSteps`` normally and ``innerSteps`` when ``dualtime == 1``.

The pure policy is:

.. math::

   \mathrm{CFL}(i) =
   \begin{cases}
     \mathrm{CFL}_{start}
       + (\mathrm{CFL}_{end}-\mathrm{CFL}_{start})\,i/n_{CFL},
       & i < n_{CFL}, \\
     \mathrm{CFL}_{end}, & i \ge n_{CFL}.
   \end{cases}

Thus the initial iteration uses ``cflst``, the ramp reaches ``cfled`` at
``ncfl``, and subsequent iterations retain the end value. The current code
does not validate a positive ramp length or constrain either CFL value.
A zero ramp length selects the end value for non-negative iteration indices.

From CFL to cell time step
--------------------------

For the CPU local-step branch, the inviscid cell estimate is

.. math::

   \Delta t_{inv,c}
     = \frac{\mathrm{CFL}\,V_c}{\lambda_{inv,c}}.

When the viscous time-step contribution is enabled, OneFLOW computes a
viscous estimate with the same cell-volume/spectral-radius structure and
combines the two estimates:

.. math::

   \Delta t_c =
     \frac{\Delta t_{inv,c}\,\Delta t_{vis,c}}
          {\Delta t_{inv,c}+\Delta t_{vis,c}}.

The spectral radii are assembled from face contributions. The viscous
contribution depends on the selected viscous spectral-radius model and is
only included when a viscous model and ``visTimeStepModel`` enable it.
The local path also applies ``max_time_ratio`` as an upper cap relative to
the minimum local step when that ratio is positive.

Time-step model dispatch
------------------------

The current unstructured NS path in ``UTimeStep::CalcTimeStep`` has three
branches:

* ``timestepModel == 0``: retain per-cell local time steps.
* ``timestepModel == 1``: fill cells with ``ctrl.pdt``, initialized from
  ``global_dt``. This branch does not derive the step from the CFL ramp.
* Other values: compute local steps, find the minimum, fill all cells with
  that minimum, and update ``ctrl.pdt``. Checked-in case comments label
  value 2 as the local-global mode; the implementation's final ``else``
  branch is broader than that documented value.

The supported HIP fast path is currently gated to local time stepping
(``timestepModel == 0``) and computes the inviscid local step. Do not
interpret it as evidence that the HIP path has the CPU viscous combination.

Source map
----------

* Parameter loading: ``Iteration::Init`` in
  ``codes/global/src/Iteration.cpp``
* Iteration selection and CFL policy call: ``TimeStep::CalcCfl`` in
  ``codes/ns/src/TimeStep.cpp``
* Pure interpolation: ``CflAtIteration`` in
  ``codes/ns/src/TimeStepPolicy.cpp``
* Time-step model dispatch and cell assembly:
  ``UTimeStep::CalcTimeStep`` in ``codes/uns/src/UTimeStep.cpp``
* Fixed/global step input: ``Ctrl::Init`` in
  ``codes/global/src/Ctrl.cpp``
* Tests: ``tests/solver/timestep_policy_test.cpp``
* Contract Workbench id: ``parameter/cfl-ramp``

What this proves / does not prove
---------------------------------

**Proves:** source-level parameter loading, iteration-counter selection,
linear interpolation, pure-policy test declarations, and the current
unstructured NS time-step dispatch and formulas.

**Does not prove:** that the linked tests have run, that every solver uses
this path, that a configured CFL is stable for a given mesh/flow, or that a
case reaches a converged or accurate solution. Units and nondimensional
interpretation of ``global_dt`` remain case/solver dependent.
