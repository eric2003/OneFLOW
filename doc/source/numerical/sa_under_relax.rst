SA under-relaxed working-variable update
========================================

One-equation SA advances the working variable with optional under-relaxation:

.. math::

   q \leftarrow (1-c)\,q_0 + c\,(q_0+\Delta q)

where \(c=1\) recovers the full explicit update.

References
----------

* Implementation practice for SA positivity / stability.
* See ``sa_update_bounds.rst`` and ``verification_layers.rst``.

Pure API: ``SaUnderRelaxedUpdate``.
