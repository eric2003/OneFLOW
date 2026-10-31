LUSGS viscous spectral contribution
===================================

After the inviscid flux Jacobian increment, the viscous spectral radius
damps the neighbor delta:

.. math::

   r \leftarrow r - \lambda_{\mathrm{vis}}\,\Delta q

Model 1 stores a per-volume radius
:math:`\lambda_{\mathrm{vis}}=\lambda_{\mathrm{face}}/V`.
Model 2 uses the face radius directly.

Face density for the viscous radius models is the arithmetic mean of the
left/right states.

References
----------

* Yoon–Jameson LUSGS; viscous spectral radius practice.
* See ``lusgs_flux_increment.rst`` and ``verification_layers.rst``.

Pure APIs: ``LusgsViscousSpectrumContribution``,
``LusgsFaceArithmeticMean``, ``LusgsViscousRadiusPerVolume``.
Wiring: ``UNsLusgs::CalcViscousTerm``.
