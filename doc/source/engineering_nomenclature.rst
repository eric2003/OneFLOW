Engineering nomenclature and semantic cross-reference
=======================================================

Purpose
-------

CFD variable names are not globally interchangeable identifiers. A useful engineering knowledge base connects names to physical concepts while retaining software, data role, units, dimensions, applicability, and source evidence.

The initial catalog cross-references CGNS/SIDS, OpenFOAM, SU2, and public ANSYS Fluent documentation. It is curated, not exhaustive, and does not infer equivalence from spelling alone.

Semantic record model
---------------------

Each concept has a stable OneFLOW ID, label, symbol, definition, dimensional exponents, reference unit, software-scoped aliases, an equivalence level (exact, conditional, or related), applicability caveats, and source IDs. Roles distinguish solution fields, mesh geometry, post-processing variables, derived fields, and input parameters. For example, SU2 FREESTREAM_PRESSURE is a freestream input parameter, not a local pressure field.

Why spelling-only matching is unsafe
-------------------------------------

OpenFOAM p is intentionally mapped to two concepts. In common incompressible solvers it denotes kinematic pressure p_k = p_s / rho, with dimensions of velocity squared. In compressible solvers p can denote dimensional static pressure. Solver and field dimensions must decide the interpretation.

Fluent Static Pressure is gauge pressure relative to operating pressure, while Absolute Pressure adds operating pressure. They are separate concepts. CGNS Pressure is the standardized static-pressure quantity, but stored values may be dimensional or normalized according to DataClass, units, and conversion metadata.

CGNS is a useful anchor because SIDS defines standardized data-name identifiers, physical meaning, and dimensional exponents. It also represents vector quantities as scalar component arrays, so VelocityX, VelocityY, and VelocityZ are components, not arbitrary aliases.

Initial covered concepts
------------------------

* Static pressure: CGNS Pressure; OpenFOAM p; Fluent Static Pressure. Check OpenFOAM scaling and Fluent gauge reference.
* Kinematic pressure: OpenFOAM p in common incompressible solvers. Dimensions m2/s2, not Pa.
* Absolute pressure: Fluent Absolute Pressure; includes operating pressure.
* Total/stagnation pressure: CGNS PressureStagnation; Fluent Total Pressure; OpenFOAM total(p). Model, frame, and derivation options matter.
* Density: CGNS Density; OpenFOAM rho; Fluent Density. SU2 FREESTREAM_DENSITY is a reference input, not a field.
* Velocity vector: CGNS VelocityX/Y/Z; OpenFOAM U; Fluent velocity components. Frame, component, and location matter.
* Mesh coordinates: CGNS CoordinateX/Y/Z; OpenFOAM points. Mesh geometry is distinct from a flow field.

* Mach number: CGNS Mach; SU2 MACH_NUMBER; OneFLOW mach_ref; Fluent Mach Number. Reference inputs and local post-processing fields are different roles.
* Reynolds number: CGNS Reynolds; SU2 REYNOLDS_NUMBER; OneFLOW reynolds. The OneFLOW read site is recorded, while its exact reference-scale convention remains a trace item.
* Static temperature: CGNS Temperature; OpenFOAM T. Preserve the distinction between thermodynamic temperature and temperature-like model inputs.
* Molecular and kinematic viscosity: CGNS ViscosityMolecular and ViscosityKinematic; OpenFOAM mu and nu where exposed. Their dimensions differ and must not be interchanged.
* Turbulence closure fields: CGNS TurbulentEnergyKinetic, TurbulentDissipation, TurbulentDissipationRate, TurbulentSANuTilde; OpenFOAM k, epsilon, omega, nut, mut, and nuTilda. Model and dimensional conventions remain part of each mapping.

Use in OneFLOW
--------------

The Qt-independent service is tools/verification/nomenclature_catalog.py, with data in tools/verification/nomenclature/cfd_core.json. Exact alias lookup returns all matches instead of silently choosing one. A future Workbench view can use the same service for parameter discovery, code navigation, formula links, and cross-solver case translation.

Expansion workflow
------------------

#. Identify a concept in equations or a OneFLOW runtime field/parameter.
#. Record OneFLOW declarations, reads, units, and storage location.
#. Add external names only when supported by a standard, source code, or manual.
#. Mark each mapping exact, conditional, or related.
#. Add tests for ambiguity and important non-equivalences.
#. Link concepts to formulas, contracts, and representative cases.

Reference sources
-----------------

* CGNS/SIDS Appendix A: https://cgns.org/standard/SIDS/appendix_a.html
* CGNS Grid Coordinates and Flow Solution: https://cgns.org/standard/SIDS/grid.html
* OpenFOAM kinematic pressure: https://doc.openfoam.com/2312/tools/processing/solvers/algorithm-kinematic-pressure/
* OpenFOAM pressure function object: https://doc.openfoam.com/2212/tools/post-processing/function-objects/field/pressure/
* SU2 configuration template: https://github.com/su2code/SU2/blob/master/config_template.cfg
* ANSYS Fluent field variable definitions: https://ansyshelp.ansys.com/public/Views/Secured/corp/v242/en/flu_ug/flu_ug_fvdefs.html
* ANSYS CFD-Post variable mapping: https://ansyshelp.ansys.com/public/Views/Secured/corp/v251/en/cfd_post/SECTION00323.html
* OpenFOAM turbulenceFields: https://doc.openfoam.com/2212/tools/post-processing/function-objects/field/turbulenceFields/
