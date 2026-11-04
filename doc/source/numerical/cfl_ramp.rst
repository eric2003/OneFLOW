Time stepping: CFL, physical time, and dual-time status
==========================================================

Purpose and evidence boundary
-----------------------------

This page is a living source-trace and verification memorandum for the
current time-stepping implementation. It records what the checked-in source
does, the mathematical relationships needed to reason about it, known
verification gaps, and the tests/cases that should be added later.

An implementation declaration is not a verification result. A linked test
declaration is not evidence that the test was built or executed. A passing
regression comparison is not, by itself, independent numerical validation.
The maturity of each claim must remain explicit.

Recommended evidence levels:

* implemented: a source path exists;
* unit-tested: an isolated policy or numerical kernel has executed tests;
* integration-tested: the connected runtime path has executed tests;
* numerically verified: results have been compared with an analytical or
  independently computed reference, with error behavior assessed;
* case-validated: a representative CFD case has been run and assessed;
* regression-tested: a previously accepted result is reproduced within a
  declared comparison policy.

These levels are complementary, not interchangeable.

CFL schedule
------------

The parameters ``cflst``, ``cfled``, and ``ncfl`` define a CFL
schedule for the unstructured Navier-Stokes time-step path. They do not
directly specify a physical time increment.

``Iteration::Init`` reads the three values from the runtime parameter
database. ``UTimeStep::CalcTimeStep`` calls ``CalcCfl()`` before selecting
the CPU or supported HIP path. The optional ``idualtime`` parameter is loaded
into ``Iteration::dualtime`` (default zero for older parameter files).
``TimeStep::CalcCfl`` selects ``outerSteps`` for physical-time stepping
and ``innerSteps`` when ``idualtime == 1``.

The pure schedule is

.. math::

   \mathrm{CFL}(i) =
   \begin{cases}
     \mathrm{CFL}_{start}
       + (\mathrm{CFL}_{end}-\mathrm{CFL}_{start})\,i/n_{CFL},
       & i < n_{CFL}, \\
     \mathrm{CFL}_{end}, & i \ge n_{CFL}.
   \end{cases}

The initial iteration uses ``cflst``, the ramp reaches ``cfled`` at
``ncfl``, and later iterations retain the end value. The current code does
not validate a positive ramp length or constrain either CFL value. A zero
ramp length selects the end value for non-negative iteration indices.

From CFL to cell time step
--------------------------

For the CPU local-step branch, the inviscid cell estimate is

.. math::

   \Delta t_{inv,c}
     = \frac{\mathrm{CFL}\,V_c}{\lambda_{inv,c}}.

When the viscous time-step contribution is enabled, the estimates are
combined as

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

Physical-time accumulation and stopping
----------------------------------------

``Ctrl::Init`` initializes ``ctrl.pdt`` from ``global_dt`` and sets
``ctrl.currTime`` to zero. It also seeds ``ctrl.pdt1 = ctrl.pdt`` for
the first history-dependent calculation.

In the audited ``MG::Run`` implementation, ``ctrl.currTime += ctrl.pdt``
is the physical-time accumulator update. Its placement depends on
``startStrategy``:

* For ``startStrategy == 2`` or ``3``, the increment occurs inside the
  inner loop before ``SolveInnerIter``; both ``outerSteps`` and
  ``innerSteps`` are incremented there.
* In the other branch, the increment occurs once per outer-loop pass before
  the inner convergence loop.

``SimuIterState::Running`` chooses the stop condition:
``outerSteps >= maxSteps`` when ``iexitflag == 0``; otherwise
``currTime >= maxTime``. Since the time increment occurs before the next
outer-loop running check, the audited path does not show a final-step clip to
``maxTime``. The actual accumulation cadence is therefore coupled to
``startStrategy`` and should not be generalized to every solver without
tracing its run loop.

This is a source observation, not a statement that either cadence is
numerically correct for every configuration. In particular, the increment
before the solve raises questions about rejected, interrupted, or
non-converged physical steps. These need explicit lifecycle tests before
changing the behavior.

Physical time, pseudo time, and history intervals
--------------------------------------------------

Use distinct notation for the physical interval and the pseudo-time
iteration:

* \(t^n\): accepted physical time at level \(n\);
* \(\Delta t_n = t^{n+1}-t^n\): current physical interval;
* \(\Delta t_{n-1}\): previous accepted physical interval;
* \(\tau\): pseudo-time used to converge the current physical level.

For a variable-step method, define the interval ratio

.. math::

   r_n = \frac{\Delta t_n}{\Delta t_{n-1}}.

For reference, the variable-step BDF2 derivative is

.. math::

   \left.\frac{dQ}{dt}\right|_{n+1}
   \approx \frac{1}{\Delta t_n}
   \left[
     \frac{1+2r_n}{1+r_n}Q^{n+1}
     -(1+r_n)Q^n
     +\frac{r_n^2}{1+r_n}Q^{n-1}
   \right].

This formula is included as a mathematical reference for history-dependent
integration. It is not a claim that every ``linearTwoStepMethods`` selector
currently implements this variable-step BDF2 formula. The exact mapping from
the legacy coefficient table to each method must be derived and independently
checked before assigning formal order or stability properties.

A generic dual-time residual equation can be written as

.. math::

   \frac{\partial Q}{\partial \tau}
   + \frac{\alpha_0 Q^{n+1}
            +\alpha_1 Q^n
            +\alpha_2 Q^{n-1}}{\Delta t_n}
   + R(Q^{n+1}) = 0.

The pseudo-time iteration solves for the current physical level. Advancing
\(\tau\) is not the same operation as accepting another physical time
level. A correct lifecycle must define when the physical solution and its
history are committed, and must not rotate physical history on every
pseudo-time or multigrid smoothing iteration.

Previous-step interval bootstrap and remaining lifecycle gap
------------------------------------------------------------

``Ctrl::Init`` seeds ``ctrl.pdt1 = ctrl.pdt`` immediately after loading
``global_dt``. This gives the first history-dependent unsteady calculation
a nonzero previous-step interval when the configured step is valid.
``TimeIntegration::CalcSrcCoeff`` uses

.. math::

   sc_1 = \frac{bsc_1}{\Delta t_n}, \qquad
   sc_3 = \frac{bsc_3}{\Delta t_{n-1}}, \qquad
   sc_2 = -sc_1-sc_3,

with the current implementation represented by ``ctrl.pdt`` and
``ctrl.pdt1``.

This is only an initialization seed. The audited source tree has no
subsequent assignment that shifts the accepted current interval into
``pdt1`` when ``pdt`` changes. In particular,
``UTimeStep::CalcLgTimeStep`` can replace ``pdt`` with the minimum local
step, while the history interval remains unchanged. The startup seed removes
the zero-initialized interval; it does not complete variable-step history
management.

Known risks and questions to retain
-----------------------------------

The following are tracked verification gaps, not all confirmed numerical
failures. They must remain visible until a test or source change resolves
them.

#. **History rotation:** there is no explicit accepted-physical-step operation
   that advances ``pdt1`` to the prior accepted ``pdt``. A changing step
   can therefore make the history coefficient use a stale interval.
#. **Time-step ownership:** local-global time-step selection writes
   ``ctrl.pdt``. It is not yet established that every consumer sees a
   consistent current interval or that history is captured at the correct
   boundary.
#. **Physical-time acceptance:** ``MG::Run`` increments ``currTime``
   before ``SolveInnerIter`` in the audited paths. Failure, early exit,
   divergence, or cancellation may leave time and solution state out of
   sync unless another layer restores or rejects the step.
#. **Dual-time evidence:** ``idualtime`` currently affects CFL iteration
   indexing, and history-dependent spectrum/source paths exist, but this is
   not evidence of an end-to-end dual-time solver lifecycle.
#. **Pseudo-time convergence:** residual thresholds, maximum pseudo-time
   iterations, multigrid cycles, and the physical-step acceptance condition
   need separate definitions and tests.
#. **Startup method:** multistep methods require a defined first-step
   bootstrap. The current initialization of ``pdt1`` alone does not prove
   that the solution history required by every method is initialized
   consistently.
#. **Invalid intervals:** divisions by ``pdt`` and ``pdt1`` require
   finite positive intervals. The relevant parameter and runtime validation
   boundaries have not yet been established.
#. **Final-time handling:** the audited loop does not show clipping the last
   physical interval to ``maxTime - currTime``. Overshoot behavior should be
   tested and specified rather than assumed.
#. **Strategy dependence:** the two ``MG::Run`` branches use different
   counter and time-accumulation cadences. Other solver paths may differ and
   require their own call-chain audit.
#. **CPU/HIP equivalence:** the supported HIP local-step path currently does
   not establish equivalence with the CPU viscous contribution.

Verification backlog
--------------------

The work should be added incrementally as real solver paths and representative
cases become available. Do not mark these items verified merely because the
contract or test declaration has been added.

Unit and policy tests:

* finite-positive interval validation, including zero, negative, NaN, and
  infinity;
* fixed-step history bootstrap and constant-step history rotation;
* variable-step history rotation and interval-ratio coefficient evaluation;
* startup behavior for each supported time-integration method;
* CFL schedule selection for physical-time and dual-time modes;
* stopping-time and final-interval clipping policy.

Runtime lifecycle and integration tests:

* one physical step with an observable accepted/rejected outcome;
* multiple physical steps with changing \(\Delta t_n\);
* pseudo-time iteration that does not advance physical time per inner
  iteration;
* multigrid smoothing/cycles that do not rotate physical history;
* early convergence, maximum-iteration exit, divergence, and cancellation;
* both audited ``startStrategy`` branches;
* local, global, and local-global time-step models;
* CPU/HIP parity for configurations claimed to be supported.

Numerical verification and case validation:

* manufactured time-dependent solution with known temporal error;
* constant-step temporal refinement and observed order;
* variable-step temporal refinement with controlled interval ratios;
* a dual-time benchmark with independently checked residual and physical
  time histories;
* at least one unsteady compressible and one unsteady incompressible
  representative case, when their runtime paths are ready;
* comparison of physical observables, not only residual convergence.

For each numerical result, retain the governing equation, initial and boundary
conditions, mesh/time-step parameters, solver configuration, reference
solution, error norm, observed convergence rate, command, revision, and
artifacts. A regression baseline should be labelled as reproducibility
evidence unless it is independently tied to a numerical reference.

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
  ``codes/global/src/Ctrl.cpp`` and the solver-specific control
  initializers
* Physical-time loop: ``MG::Run`` in
  ``codes/multigrid/src/Multigrid.cpp``
* History-dependent coefficients: ``TimeIntegration::CalcSrcCoeff`` in
  ``codes/unsteady/src/TimeIntegration.cpp``
* Dual-time spectrum consumers: ``UNsSpectrum.cpp``,
  ``UINsSpectrum.cpp``, and ``UTurbSpectrum.cpp``
* Existing CFL policy tests: ``tests/solver/timestep_policy_test.cpp``
* Existing parameter contract: ``parameter/cfl-ramp``

What this proves / does not prove
---------------------------------

**Proves:** source-level parameter loading, iteration-counter selection,
linear interpolation, dual-time iteration-index selection, pure-policy test
declarations, current unstructured NS time-step dispatch, and the formulas
documented above as implementation/reference relationships.

**Does not prove:** that every linked test has run; that every solver uses
this path; that variable-step history is rotated correctly; that dual-time
integration is exercised end-to-end; that a configured CFL is stable for a
given mesh/flow; that the selected time integrator has its claimed formal
order; or that a case reaches a converged or accurate solution. Units and
nondimensional interpretation of ``global_dt`` remain case/solver
dependent.

Evolution policy
----------------

This memorandum is intentionally allowed to contain unresolved items.
As refactoring and case coverage expand, each item should move through
explicit states such as ``observed``, ``test_planned``,
``test_declared``, ``executed``, ``numerically_verified``,
and ``resolved``. Keep the original observation and link the resolving
commit and evidence; do not erase the history when implementation changes.
