Boundary Condition Architecture
=================================

Purpose
-------

Boundary conditions are a first-class physical and numerical domain in
OneFLOW. They connect the physical problem to the mathematical model, mesh
boundary entities, discretization, solver implementation and runtime state
data. They are not merely mesh labels or converter string substitutions.

This page records current Gmsh-to-CGNS-to-OneFLOW type mappings and the
information required for a richer, loss-aware semantic contract. A mapped
type is not proof that a compatible CFD model or boundary state is configured.

Three layers that must remain distinct
--------------------------------------

1. Boundary identity identifies a patch or interface in a mesh: Inlet,
   Outlet, UpperWall, Cylinder, PeriodicA or Fluid. Names, source IDs,
   dimensions, entity membership, orientation and provenance belong to Mesh.
2. Physical boundary semantics states what the boundary means in the
   mathematical problem: prescribed velocity, mass flow, static pressure,
   characteristic inflow, pressure outlet, no-slip wall, inviscid slip wall,
   isothermal wall, specified heat flux, symmetry, farfield, overset or
   periodic interface.
3. Numerical realization and state data describe how a solver enforces
   the condition: ghost-cell reconstruction, flux constraints, characteristic
   relations, viscous wall treatment, thermal closure, turbulence variables,
   reference states, gradients and normal direction.

.. code-block:: text

   Mesh patch identity
       (name, source ID, dimension, members, provenance)
                         |
                         v
   Physical boundary semantic contract
       (condition family, mathematical constraints, required state)
                         |
              +----------+----------+
              |                     |
              v                     v
       Format representation   Solver realization
       Gmsh Physical Group     OneFLOW BC category,
       CGNS BC_t / BCType_t    numerical treatment and runtime data
              |                     |
              +----------+----------+
                         v
                Verification evidence

The semantic contract is the bridge. A format adapter must not invent
physical data, and a solver adapter must not claim a physical condition was
realized merely because it received a broad boundary ID.

Common physical condition families
-----------------------------------

* Inlet: velocity, mass-flow, total-condition, pressure or characteristic
  inflow. The formulation determines which state variables are prescribed.
* Outlet: static-pressure, extrapolative, characteristic or other outflow.
  Subsonic and supersonic formulations have different data needs.
* Wall: impermeability and tangential velocity constraints, distinguishing
  inviscid slip, viscous no-slip, isothermal and specified heat-flux treatments.
* Symmetry: zero normal velocity and appropriate normal-gradient/flux
  constraints for the governing equations.
* Farfield: free-stream/reference state and a formulation appropriate to
  local characteristic direction and flow regime.
* Interface: internal connectivity, periodic transforms, overset
  interpolation or another coupling contract.
* Special geometric boundaries: axis, pole, degenerate edge/point and wedge
  constraints, dependent on dimension and coordinate system.

These are semantic families, not a claim that every formulation is
implemented by every OneFLOW solver.

Gmsh Physical Group name to CGNS type
-------------------------------------

The current adapter removes non-alphanumeric characters and compares
case-insensitively. The original patch name remains a separate boundary
identity; the mapping selects a coarse CGNS type.

.. list-table::
   :header-rows: 1
   :widths: 32 24 44

   * - Gmsh Physical Group aliases
     - CGNS BCType_t
     - Meaning and limitation
   * - Inlet, Inflow, VelocityInlet, MassFlowInlet, PressureInlet
     - BCInflow (9)
     - Generic inflow. Alias spelling does not encode prescribed state variables.
   * - Outlet, Outflow, PressureOutlet
     - BCOutflow (13)
     - Generic outflow. Target pressure and characteristic treatment are separate.
   * - Symmetry, SymmetryPlane
     - BCSymmetryPlane (16)
     - Symmetry-plane representation.
   * - Farfield, Freestream
     - BCFarfield (7)
     - Farfield representation; free-stream state is separate.
   * - NoSlipWall, ViscousWall
     - BCWallViscous (22)
     - Viscous-wall type; thermal treatment is not specified.
   * - SlipWall, InviscidWall
     - BCWallInviscid (21)
     - Inviscid-wall type.
   * - Wall, LowerWall, UpperWall, BodyWall
     - BCWall (20)
     - Broad wall type; do not infer slip/no-slip or thermal behavior.
   * - Unrecognized names
     - BCTypeUserDefined (1)
     - Unknown semantics remain unknown; do not guess.

Physical Group names are user-defined, not a standardized Gmsh physical
boundary vocabulary. PressureInlet and VelocityInlet both map to BCInflow.
Their distinction must remain in semantic metadata and boundary-state data.
Likewise, LowerWall and UpperWall identify different patches but do not state
their physical wall treatment.

CGNS BCType_t reference
----------------------

The following numeric values are the CGNS header enumeration values, not
OneFLOW boundary IDs.

.. list-table::
   :header-rows: 1
   :widths: 8 36 56

   * - Value
     - CGNS enum
     - Type-level meaning
   * - 0
     - BCTypeNull
     - Null / no boundary classification.
   * - 1
     - BCTypeUserDefined
     - User-defined classification.
   * - 2
     - BCAxisymmetricWedge
     - Axisymmetric wedge.
   * - 3
     - BCDegenerateLine
     - Degenerate line.
   * - 4
     - BCDegeneratePoint
     - Degenerate point.
   * - 5
     - BCDirichlet
     - Dirichlet-type condition.
   * - 6
     - BCExtrapolate
     - Extrapolation.
   * - 7
     - BCFarfield
     - Farfield.
   * - 8
     - BCGeneral
     - General condition.
   * - 9
     - BCInflow
     - Inflow.
   * - 10
     - BCInflowSubsonic
     - Subsonic inflow.
   * - 11
     - BCInflowSupersonic
     - Supersonic inflow.
   * - 12
     - BCNeumann
     - Neumann-type condition.
   * - 13
     - BCOutflow
     - Outflow.
   * - 14
     - BCOutflowSubsonic
     - Subsonic outflow.
   * - 15
     - BCOutflowSupersonic
     - Supersonic outflow.
   * - 16
     - BCSymmetryPlane
     - Symmetry plane.
   * - 17
     - BCSymmetryPolar
     - Polar symmetry.
   * - 18
     - BCTunnelInflow
     - Tunnel inflow.
   * - 19
     - BCTunnelOutflow
     - Tunnel outflow.
   * - 20
     - BCWall
     - Wall.
   * - 21
     - BCWallInviscid
     - Inviscid wall.
   * - 22
     - BCWallViscous
     - Viscous wall.
   * - 23
     - BCWallViscousHeatFlux
     - Viscous wall with prescribed heat flux.
   * - 24
     - BCWallViscousIsothermal
     - Isothermal viscous wall.
   * - 25
     - FamilySpecified
     - Family-specified condition.

A CGNS BC node also has a grid location and point/element list, and may carry
additional BC data or datasets. Correctly writing BCType_t alone does not
establish that patch entities, normals, location, state arrays or units are
correct. Periodic interfaces are normally represented through zone
connectivity and periodic transform metadata, not inferred from a name alias
or a BCType_t member.

CGNS type to OneFLOW category
-----------------------------

The current BcTypeMap::Init implementation explicitly maps:

.. list-table::
   :header-rows: 1
   :widths: 48 28 24

   * - CGNS type(s)
     - OneFLOW category
     - Current meaning
   * - BCTypeNull
     - BC::INTERFACE
     - Interface classification.
   * - BCTypeUserDefined
     - BC::GENERIC_2
     - Generic fallback.
   * - BCSymmetryPlane
     - BC::SYMMETRY
     - Symmetry.
   * - BCInflow, BCInflowSubsonic, BCInflowSupersonic, BCTunnelInflow
     - BC::INFLOW
     - Inflow variants collapse to one ID.
   * - BCOutflow, BCOutflowSubsonic, BCOutflowSupersonic, BCTunnelOutflow
     - BC::OUTFLOW
     - Outflow variants collapse to one ID.
   * - BCWall, BCWallInviscid, BCWallViscous, BCWallViscousHeatFlux, BCWallViscousIsothermal
     - BC::SOLID_SURFACE
     - Wall subtypes collapse to one broad ID.
   * - BCDegenerateLine
     - BC::POLE
     - Pole / degenerate line.
   * - BCFarfield
     - BC::FARFIELD
     - Farfield.
   * - BCGeneral
     - BC::OVERSET
     - General / overset.

The reverse mapping is not perfectly symmetric. BC::EXTRAPOLATION maps to
BCExtrapolate on export, but BCExtrapolate is not explicitly mapped on import
and reaches the generic fallback. Other CGNS types not listed above also use
the current default BC::GENERIC_2. This is a coverage gap, not a semantic
equivalence.

Information that must be preserved
----------------------------------

A robust boundary contract must retain, where applicable:

* stable patch identity, original name, source format, source entity IDs and
  dimension;
* physical condition family and formulation (pressure inlet versus velocity
  inlet, for example);
* governing-equation/model applicability and flow-regime assumptions;
* prescribed values, functions/profiles, units, reference quantities and
  coordinate frame;
* required and derived variables, including thermodynamic and turbulence
  quantities where relevant;
* wall thermal model (adiabatic, isothermal, heat flux or another declared
  treatment), wall motion and roughness where supported;
* face/edge membership, grid location, orientation, outward normal and
  interface/periodic connectivity;
* mapping decisions, source/target enum values, adapter version and explicit
  loss/unsupported diagnostics;
* solver implementation, parameter consumers, run configuration and
  verification artifacts.

Missing data is not a default value. Unknown semantics remain unknown until
supplied or deliberately selected.

Contract ownership and verification lifecycle
----------------------------------------------

Boundary Conditions owns physical semantics and format/solver mapping
knowledge. Mesh owns geometric patch identity, entity membership, topology
and orientation. Parameter owns runtime settings and consumers. Case owns
the physical scenario and selected conditions. Runtime owns when and how the
solver consumes them. Workflow connects conversion and solver steps through
explicit artifact contracts.

The intended evidence chain is:

1. declare a physical boundary contract and required state data;
2. bind it to mesh patch identities without changing those identities;
3. map it to source/target format representations and report information loss;
4. resolve it to a supported OneFLOW category and solver treatment;
5. validate required state data before launching the solver;
6. run a representative numerical/physical case and retain logs, outputs,
   comparisons and revision information.

A successful format conversion verifies only the conversion claims it checks.
It does not prove solver compatibility, correct physical state, numerical
stability or physical accuracy. Mesh-to-native-grid validation, boundary
semantic validation and CFD physical validation are separate evidence levels.

Workbench presentation
----------------------

Boundary Conditions is a first-class entry on the OneFLOW Workbench home
screen and a first-class domain in Contract Workbench. Its shared navigation
covers Overview, Principle, Reference, Architecture, Implementation, Code,
Documentation and Verification / Evidence. Its domain-specific Boundary Model
view presents semantic and format/solver mapping contracts. It is not
embedded in the Mesh Validation dialog.

Mesh continues to describe boundary entities and topology. Mesh Validation
invokes conversion and verification services and links to Boundary Conditions
contracts; it does not own the semantic catalog.

References
----------

* CGNS SIDS: https://cgns.github.io/CGNS_docs_current/sids/
* CGNS Mid-Level Library: https://cgns.github.io/CGNS_docs_current/midlevel/
* Gmsh Reference Manual: https://gmsh.info/doc/texinfo/gmsh.html
* Implementation: codes/geometry/src/Boundary.cpp,
  codes/geometry/src/GridConversion.cpp
* mesh_architecture
* mesh_to_cfd_validation
* contract_workbench_information_architecture


Dedicated Boundary Conditions Workbench
---------------------------------------

The home-screen Boundary Conditions entry opens a dedicated domain workbench. Contract Workbench retains a separate Boundary Conditions contract domain for cross-domain governance; it is not the primary knowledge browser.

The versioned data catalog is tools/verification/boundary_conditions/registry.json. Stable semantic IDs are independent of CFL3D type numbers and OneFLOW enum values. Entries describe physical formulation, equations, parameters, data location, format mapping, implementation status and verification state.

The first curated entries cover CFL3D Version 5.0 types 1002, 1003, 2002, 2003 and 2004. The manual version and section are retained with each entry. Compact equations in the catalog are reference summaries; complete coordinate conventions, assumptions, normalization and derivations must be transcribed before implementing an equivalent solver treatment.

Mappings are qualified, not treated as equivalences. CGNS BCOutflow and OneFLOW OUTFLOW do not encode a prescribed pressure ratio. SOLID_SURFACE does not distinguish wall velocity or thermal closure. Such information must remain in the semantic definition and parameter payload.

The GUI is read-only in this first slice. It does not configure solver runs or claim numerical support. not_verified and category_only are explicit capability states. Future stages include exact equation transcription, typed parameter validation, solver capability registration and executable formula/numerical verification.
