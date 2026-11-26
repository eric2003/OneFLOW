Turbulence far-field state selection
====================================

Turbulence variables at far-field faces follow the mean-flow normal
direction: interior values when fluid leaves the domain, freestream
turbulence when fluid enters.

.. math::

   \text{use interior} \iff v_n \ge 0

Supersonic: \(v_n\) is the interior normal velocity. Subsonic: \(v_n\) is
the Riemann-combined boundary normal velocity (same acoustic practice as
the NS far-field).

References
----------

* Characteristic far-field practice; NASA TMR turbulence BC notes.
* See ``turb_wall_bc.rst`` and ``farfield_riemann.rst``.

Pure API: ``TurbFarfieldUsesInterior``.
