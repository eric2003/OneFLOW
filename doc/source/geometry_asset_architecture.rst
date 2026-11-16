Geometry Asset Architecture
============================

Purpose
-------

A geometry asset is not a file extension and is not synonymous with STL. It is
the reusable identity and design intent of a geometric object, together with
its definitions, provenance and available representations.

One asset may have several representations. For example, a sphere may be
defined by an analytic equation, instantiated by a parameter set, represented
as an exact CAD surface, tessellated into STL, and discretized into a volume
mesh. These are related artifacts, not unrelated geometry records.

Asset and representation
------------------------

The asset definition describes what the geometry is. A representation
describes how a particular tool or workflow stores or approximates it.

The initial representation roles include:

* native definition: equations, parameters or a source model;
* CAD B-rep or CAD surface: exact curves, surfaces, topology and solids;
* surface mesh: polygonal or triangular boundary approximation;
* volume mesh: discretized computational domain;
* point cloud: sampled geometric observations;
* implicit field: signed-distance or other scalar/vector field definition;
* visualization and exchange: renderable or interoperable derivatives.

Examples of formats include STEP/IGES/BREP for CAD exchange or kernel models;
STL/OBJ/PLY/3MF for tessellated surfaces; MSH/CGNS/VTU for meshes; and
application-native or scientific-data formats for other representations.
These are format categories, not a claim that every importer is implemented.
STEP and IGES carry richer CAD boundary information than STL's triangle
surface representation. A mesh-to-CAD conversion must not be described as
reconstructing the original exact design intent unless that reconstruction is
actually performed.

Definition kinds
----------------

The registry distinguishes analytic, parametric, procedural, constructive,
imported, measured and implicit definitions. These kinds describe how the
asset's geometry is defined, not which file extension is used.

A parameterized cylinder, for example, remains the same geometry asset when
its surface tessellation density changes. The changed STL is a new
representation generated with different discretization parameters; the
analytic radius and length are not changed.

Identity and provenance
-----------------------

Every asset has a stable ID independent of its filenames and output folders.
Its provenance should link to the strongest available source:

* mathematical equations and parameter definitions;
* tutorial or development Case that demonstrates the construction;
* generator source path, algorithm and version;
* external source or user-supplied origin where applicable;
* source revision and content checksum when known.

A derived representation records its format, path, generator, status,
checksum when available, and parent representation. This supports a lineage
such as:

.. code-block:: text

   Geometry Asset
     -> analytic / parametric definition
     -> CAD or procedural source
     -> surface representation (e.g. STL)
     -> volume mesh (e.g. MSH)
     -> converted mesh (e.g. CGNS)
     -> solver input and validation evidence

Not every asset follows every stage. A measured point cloud may be the primary
source; a CAD asset may go directly to a mesher; a mesh may be imported without
an available exact CAD source. Missing provenance must be stated as unknown,
not filled with an assumed origin.

Initial repository registry
---------------------------

The registry is located at `grid/geometry/`:

* `asset.schema.json` defines the versioned asset metadata contract.
* `catalog.json` registers the initial unit sphere, 6:1 prolate spheroid and
  cylinder internal-flow domain.
* `README.md` explains registry ownership and boundaries.

The initial entries refer to the existing deterministic generator
`grid/stl_to_gmsh/generate_stl_geometry.py` and the geometry construction
documentation in `stl_geometry_generation.rst`. This is provenance for the
current reproducible generation path. It is explicitly not a claim that
historical checked-in STL fixtures were emitted by that generator.

The existing `grid/case_catalog.json` remains the runnable Case catalog.
The registry does not replace it. Likewise, `grid/stl_to_gmsh/` remains a
workflow and example collection rather than the canonical asset registry.

Format and capability policy
----------------------------

The registry can describe formats beyond those currently supported by
OneFLOW. Import, export, validation and conversion support must be reported by
the corresponding adapter. The registry itself must not load a CAD kernel or
assume that a file extension guarantees a valid or exact representation.

The first implementation is metadata-first. The deterministic STL
generator writes a companion ``.geometry.json`` manifest containing the stable
asset ID, a snapshot of the analytic/parametric definition, definition SHA-256,
generation parameters, generator name/version/source path/source SHA-256, and
the generated STL SHA-256. Definition parameters are kept separate from
representation-generation parameters: changing tessellation density creates a
new surface representation, not a new geometric identity. The legacy
``parameters`` field remains for compatibility; new consumers should use
``definition.parameters`` and ``generation_parameters``.

The source checksum identifies the exact generator source bytes used for the
run; it is not a Git revision. A Git revision is recorded only when the caller
can supply a verified revision. Historical assets must not inherit this
manifest retroactively: their original generator and revision remain unknown
unless independent evidence establishes them.

The remaining integration work is to let Geometry Library launch the
asset-specific parameter editor, carry the manifest through STL-to-Gmsh and
MSH-to-CGNS, and display the derivation chain in Mesh-to-CFD Validation.

References
----------

* :doc:`stl_geometry_generation` — current parameterized surface generator,
  mathematical definitions, patch identity and generated-artifact policy.
* :doc:`mesh_architecture` — mesh representations, snapshots and adapters.
