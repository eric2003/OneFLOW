Boundary velocity images
========================

Symmetry and viscous-wall boundary conditions build ghost velocities from
simple linear images of the interior state.

Symmetry (reflection)
---------------------

With unit face normal :math:`\mathbf{n}` and face speed :math:`v_{fn}`:

.. math::

   v_n^{\mathrm{rel}} = \mathbf{n}\cdot\mathbf{v} - v_{fn}

.. math::

   \mathbf{v}'
   = \mathbf{v} - 2\,\mathbf{n}\,v_n^{\mathrm{rel}}

No-slip wall image
------------------

.. math::

   \mathbf{v}'
   = -\mathbf{v} + 2\,\mathbf{v}_{\mathrm{wall}}

Face primitive average for the symmetry face state is the arithmetic mean
of interior and reflected states.

References
----------

* Blazek, J. *Computational Fluid Dynamics: Principles and Applications*
  (symmetry / wall boundary treatment).
* See ``references.rst`` for the citation policy.

Pure APIs: ``NsBcPolicy.h``. Tests: ``bc_policy_test``.
