Verification and Evidence Lifecycle
====================================

Purpose
-------

OneFLOW is a scientific-computing system. A successful process exit is
necessary for many checks, but it is not sufficient evidence that a numerical
method is correct or that a simulation represents physical reality. This
guide defines the evidence chain used by the Contract Workbench and explains
how it connects software engineering practice with scientific verification.

The lifecycle is intentionally layered:

.. code-block:: text

   physical model and assumptions
              |
              v
   mathematical formulation
              |
              v
   discrete method and boundary treatment
              |
              v
   implementation contract
              |
              v
   declared source / test / case / reference
              |
              v
   static source / CMake declaration resolution
              |
              v
   configured-build CTest discovery
              |
              v
   executed run + revision + configuration
              |
              v
   numerical comparison and uncertainty assessment
              |
              v
   reviewed claim status and updated documentation

Each arrow is a relationship that needs evidence. A declaration is not an
execution. A source declaration is not a built test. A successful test is not
automatically proof of every claim in its contract. A regression match is not
independent physical validation.

Evidence states and their meaning
---------------------------------

.. list-table::
   :header-rows: 1
   :widths: 20 35 45

   * - Layer
     - What it establishes
     - What it does not establish
   * - Contract-declared
     - A contract names a source, test, case, or reference.
     - That the path exists or the item was executed.
   * - Resolved
     - A supported static resolver found the declared artifact or source declaration.
     - Build success, test discovery, or runtime success.
   * - CMake-declared
     - Source-level CMake parsing finds a target and a gtest_discover_tests() declaration.
     - That CMake executed the declaration or that CTest lists the test.
   * - Discovered
     - The configured build tree's CTest JSON lists the test for the selected configuration.
     - That the test was executed or passed.
   * - Executed
     - A run record identifies command, revision, platform, outcome, and artifacts.
     - That the numerical result is scientifically accurate.
   * - Passed
     - The declared test or comparison satisfied its encoded acceptance criteria.
     - That all contract claims, configurations, or physical regimes are verified.
   * - Fresh
     - Evidence was produced from the revision being assessed.
     - That the test scope is adequate or the acceptance criteria are valid.

CTest's ``--show-only=json-v1`` is the configured-build discovery source;
it lists tests without executing them and provides command/configuration and
backtrace metadata. The source-level CMake parser is only a navigation hint.
The Workbench keeps these layers separate. A passing test record must not
promote an entire contract to verified unless the contract's full claim set
and verification obligations have been evaluated.

Verification, validation, and regression
-----------------------------------------

Use the terms precisely:

* **Code verification** asks whether the implementation correctly represents
  the intended mathematical or conceptual model.
* **Solution verification** estimates numerical error in a particular
  calculation, including iterative, spatial, and temporal discretization error.
* **Validation** assesses whether the model represents physical observations
  for a stated application domain, accounting for experimental uncertainty.
* **Regression testing** checks whether a result remains consistent with an
  accepted baseline under a controlled configuration. It is valuable, but a
  baseline match alone is neither independent verification nor validation.

This distinction is central to interpreting OneFLOW case results. A solver
case can pass output comparison while still lacking grid-convergence evidence,
independent analytical comparison, or experimental validation.

Numerical evidence should grow in strength
------------------------------------------

For an algorithm, begin with the narrowest meaningful checks and increase
scientific scope deliberately:

#. Algebraic identities, limiting cases, and local unit tests.
#. Component tests for reconstruction, fluxes, gradients, source terms, and
   boundary treatment.
#. Manufactured solutions or exact analytical solutions where available.
#. Conservation and consistency checks.
#. Iterative convergence and residual behavior.
#. Spatial and temporal refinement studies.
#. Comparison with trusted benchmark data and experimental measurements.

For a family of grids with refinement ratio :math:`r`, an observed order may
be estimated from errors :math:`E_h` and :math:`E_{h/r}`:

.. math::

   p_{obs} = \frac{\ln(E_h/E_{h/r})}{\ln(r)}.

This two-grid expression is only a diagnostic under restrictive assumptions.
A defensible observed-order and discretization-uncertainty assessment normally
requires at least three systematically refined grids, consistent quantities
of interest, converged iterative solutions, and evidence that the results are
in the asymptotic range. Do not report an order estimate without stating the
refinement ratio, error norm, grid family, solver tolerances, and relevant
assumptions.

A verification-curve page should identify the governing equation, reference
solution, axes and units, discretization parameters, error definition, data
provenance, and what the comparison does not prove.

Reference policy
----------------

Cite external theory when it materially supports an equation, numerical
method, uncertainty procedure, physical model, or verification claim. A
reference should be specific enough that another developer can locate the
original source.

Prefer, in order:

#. Original peer-reviewed papers or standards for the method or procedure.
#. Authoritative technical manuals from recognized scientific organizations.
#. Textbooks for established derivations and broader context.
#. Maintained software manuals for implementation-specific behavior.
#. Secondary tutorials as an accessible introduction, not the sole authority
   for a high-impact scientific claim.

For each important formula, record the source and distinguish it from
OneFLOW's own choices: notation, assumptions, boundary conditions, limiting
cases, discretization, stabilization, and implementation details. A citation
does not transfer correctness to the implementation; the implementation still
needs tests and numerical evidence. Respect licenses and do not copy large
sections of external manuals.

Recommended starting references
-------------------------------

* NASA NPARC Alliance, Tutorial on CFD Verification and Validation:
  https://www.grc.nasa.gov/www/wind/valid/tutorial/tutorial.html
  Useful for code verification, solution verification, validation, iterative
  convergence, consistency, and grid convergence.
* Ismail B. Celik et al., "Procedure for Estimation and Reporting of
  Uncertainty Due to Discretization in CFD Applications", Journal of Fluids
  Engineering, 130(7), 078001, 2008.
  DOI: https://doi.org/10.1115/1.2960953
  A reference procedure for reporting discretization uncertainty.
* NASA, Verification Assessment:
  https://www.grc.nasa.gov/www/wind/valid/tutorial/verassess.html
  and Spatial Grid Convergence:
  https://www.grc.nasa.gov/www/wind/valid/tutorial/spatconv.html
  Practical explanations of code and solution verification.
* Breathe, Sphinx and Doxygen integration:
  https://github.com/breathe-doc/breathe
  A possible bridge for publishing selected C++ API comments in Sphinx.

The reference list is a starting point, not an endorsement of every procedure
for every OneFLOW model. The applicable physical model and numerical method
must determine the relevant literature.

Source comments and generated API documentation
------------------------------------------------

Use source comments for information closest to the code:

* public API intent and invariants;
* units, ownership, lifetime, and mutation rules;
* preconditions, postconditions, side effects, and failure behavior;
* non-obvious numerical derivation steps and references to the governing
  equation or method.

Do not duplicate the full architecture or mathematical derivation in a header.
A useful long-term pipeline is:

.. code-block:: text

   C++ public headers
          |
          v
   Doxygen XML (selected documented API)
          |
          v
   Breathe bridge
          |
          v
   Sphinx API reference
          |
          +------> curated architecture and theory guides
          |
          +------> contract catalog and verification workbench

Sphinx is already used for curated documentation under doc/source. Doxygen
and Breathe should be introduced incrementally: first document a small,
stable public API slice, inspect the generated output, then expand. Do not
enable project-wide extraction until comment quality, exclusions, generated
file handling, and CI build cost are understood. Generated API pages
complement rather than replace reviewed architecture, algorithm, and user
documentation.

Diagrams and figures
--------------------

Prefer a diagram when it explains a relationship that is difficult to infer
from prose:

* lifecycle/state diagrams for ownership and cleanup;
* dataflow diagrams for configuration, mesh, and solver stages;
* topology sketches for entity identity and orientation;
* annotated plots for convergence, residuals, and analytical comparisons.

Every figure should have a caption, define symbols or colors, and identify
whether it describes current behavior or a proposed design. Keep source files
for generated figures and plots beside the relevant documentation or in a
documented generator. Avoid screenshots as the only record of a technical
contract.

Contract Workbench navigation
------------------------------

The Workbench is intended to make the evidence chain explorable:

* Contract to implementation and documentation.
* Contract to declared test source.
* Test to the existing typed runner.
* Run to immutable command, revision, configuration, outcome, and artifacts.
* Numerical result to formula, reference curve, error metric, and limitations.

The current Open declared test sources action is a navigation step only.
Execution remains owned by the existing test runner. A future Run action must
reuse that runner, record the exact invocation, and preserve the distinction
between a test passing and a scientific claim being established.

Evolution review
----------------

At each substantial change, record:

* the system-level risk reduced and why this area outranks alternatives;
* the code and architecture boundary changed;
* the exact contract claims affected;
* the tests, CI jobs, and numerical experiments actually executed;
* revision, compiler/backend/configuration, and result artifacts;
* limitations that remain and the next highest-value uncertainty.

The intended outcome is a system in which documentation, code, tests, and
scientific evidence continuously correct one another, rather than four
separate collections maintained by hand.
