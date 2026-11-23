SA production strain measures
=============================

Spalart–Allmaras production uses a strain/vorticity magnitude. Variants
replace pure vorticity with a strain-based or rotation-corrected measure.

Classical
---------

.. math::

   S = \sqrt{\lvert 2\,\mathrm{work}_2 + \mathrm{work}_3
   - \tfrac23 \mathrm{work}_1^2 \rvert}

Rotational corrections (implementation options 4 and 5) reduce production
in laminar vortex regions; see pure APIs for the exact algebraic form.

References
----------

* Spalart–Allmaras and rotational production corrections literature.
* NASA TMR SA.
* See ``sa_source.rst`` and ``verification_layers.rst``.

Pure APIs: ``SaStrainProductionClassical``,
``SaStrainProductionRotationalV4``, ``SaStrainProductionRotationalV5``.
