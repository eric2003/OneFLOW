OneFLOW Time-Integration Control
================================

Scope
-----

This page records the current implementation contract for the legacy
``time_integral`` and ``rk_stage`` parameters. It is a source-trace document,
not a claim that the complete temporal discretization has been formally
verified.

Selector values and dispatch
----------------------------

``codes/solver/include/TimeIntegral.h`` declares:

* ``MULTI_STAGE = 1``;
* ``LUSGS = 2``;
* ``SIMPLE = 3``.

``TimeIntegral::Init`` in ``codes/solver/src/TimeIntegral.cpp`` selects the
``RungeKutta`` function for ``MULTI_STAGE``, the ``Simple`` function for
``SIMPLE``, and the ``Lusgs`` function for the remaining values. The current
function does not explicitly reject an unknown selector value; the final
branch routes it to LUSGS. This is observed behavior, not a validation
recommendation.

The context-aware ``TimeIntegral::Relaxation`` dispatches each cycle using the
same three-way selection. The non-context overload initializes the function
pointer and invokes it for each cycle.

Stage count and coefficient schedule
------------------------------------

The global, compressible-NS, incompressible-NS, and turbulence control
initializers each read ``rk_stage`` and construct:

.. math::

   a_s = \frac{1}{s+1}, \qquad s=0,\ldots,N-1,
   \qquad N=\texttt{rk\_stage}

The implementation stores these values in ``rk_coef``. On the finest grid,
``TimeIntegral::RungeKutta`` loops over that vector and assigns
``ctrl.lhscoef = ctrl.rk_coef[stage]`` before running the residual, update,
boundary, and (for multi-zone runs) interface-exchange tasks. On a coarser
grid, the current implementation uses one legacy task sequence with
``ctrl.lhscoef = 1.0``.

These values are the coefficient schedule currently consumed by this code
path. They must not be mistaken for a complete Butcher tableau or used alone
to infer formal temporal order or stability.

Parameter vocabulary and knowledge mapping
----------------------------------------

The configuration key ``time_integral`` is a discrete method selector. Its
common descriptive aliases are "time-integration method", "time integrator",
and "temporal scheme selector". These are vocabulary aliases, not separate
OneFLOW keys. The selector values are declared in ``TimeIntegral.h``:
``MULTI_STAGE = 1``, ``LUSGS = 2``, and ``SIMPLE = 3``. The dispatch chooses
RungeKutta for ``MULTI_STAGE``, ``Simple`` for ``SIMPLE``, and Lusgs for every
other value. The OneFLOW ``SIMPLE`` label names ``TimeIntegral::Simple``; it
must not be conflated with the SIMPLE pressure-velocity coupling algorithm.

The key ``rk_stage`` is commonly described as "number of Runge-Kutta stages"
or "stage count". The coefficient schedule is ``a_s = 1/(s+1)`` for
zero-based stage indices. These descriptive names are not asserted to be
CGNS-standard parameter names: CGNS defines CFD data-model and field
nomenclature, and no direct standardized equivalent for these OneFLOW
control keys is claimed here.

The Parameter Workbench exposes these as separate records: ``time_integral``
selects the algorithm path, while ``rk_stage`` sizes the legacy coefficient
schedule used by that path. Their shared Formula & Knowledge reference
``time_integral_rk_stage_weight`` is a navigation and theory link, not proof
of a complete Butcher tableau or temporal order.

Resident Euler backend eligibility
----------------------------------

The context-aware Runge-Kutta path calls
``EvaluateEulerRungeKuttaCapability`` before using the resident Euler backend.
The current capability decision requires all of the following:

* the compressible Navier-Stokes solver;
* one local zone and the finest, single-grid level;
* five equations and the configured Lax-Friedrichs inviscid scheme;
* explicit multi-stage time integration;
* no viscous terms, source terms, or limiter;
* no interface exchange;
* a backend state that supports the requested advance.

When the decision rejects the request, or the stage coefficient vector is
empty, the implementation falls back to the legacy Runge-Kutta task sequence.
This is an acceleration-path capability gate, not a general rejection of
Runge-Kutta for other solver configurations.

Evidence and limitations
------------------------

The source-level evidence is:

* ``Ctrl.cpp``, ``NsCtrl.cpp``, ``INsCtrl.cpp``, and ``TurbCtrl.cpp`` for
  parameter reads and coefficient construction;
* ``TimeIntegral.cpp`` for selector dispatch, stage execution, and fallback;
* ``EulerRungeKuttaCapability.h/.cpp`` for the specialized resident-backend gate.

The Formula & Knowledge entry ``time_integral_rk_stage_weight`` links these
sources and this page. Its formula view deliberately has no theory curve.
Current source evidence proves the parameter-to-stage-weight wiring and the
branch structure, but does not prove a complete Butcher tableau, formal order,
stability region, convergence rate, or CFD accuracy. Those claims need
dedicated numerical references, analytical tests, and representative case
verification.
