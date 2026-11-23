Relative face-normal velocity
=============================

Purpose
-------

Evaluate the relative face-normal velocity used by inviscid flux schemes,
including moving-mesh / ALE terms.

Formula
-------

.. math::

   v_n = \mathbf{n}\cdot\mathbf{u} - v_g

Source map
----------

* ``RelativeFaceNormalVelocity`` in ``NsInvFluxPolicy``
* Downstream: ``EulerConvectiveFlux``, ``FaceMachNumber``, Roe eigenvalues
* Tests: ``tests/solver/inv_flux_policy_test.cpp``
* Theory table: ``tools/verification/data/relative_face_normal_velocity_theory.csv``
* Contract Workbench id: ``relative_face_normal_velocity``

.. image:: figures/relative_face_normal_velocity_theory.png
   :alt: Relative face-normal velocity
   :width: 520

What this proves / does not prove
---------------------------------

**Proves:** relative normal-velocity algebra matches the pure API.

**Does not prove:** full flux residual or mesh-motion coupling on a case.

References
----------

* Standard finite-volume face kinematics; ALE relative velocity.
* Related: :doc:`face_mach_number` (catalog), Euler convective flux.
