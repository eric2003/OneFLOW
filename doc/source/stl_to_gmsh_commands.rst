STL to Gmsh: 3-D Volume Mesh Command
======================================

Generate a tetrahedral volume mesh directly from one closed STL surface.
No OneFLOW executable or Case staging is needed; Python and Gmsh are required.

One mesh size
-------------

.. code-block:: console

   python grid/stl_to_gmsh/generate_volume_mesh.py --stl path/to/sphere.stl --output "$HOME/.oneflow/mesh-generation/sphere.msh" --mesh-size 0.20

Mesh-size range
---------------

.. code-block:: console

   python grid/stl_to_gmsh/generate_volume_mesh.py --stl path/to/sphere.stl --output "$HOME/.oneflow/mesh-generation/sphere.msh" --mesh-size-min 0.10 --mesh-size-max 0.20

The first command sets both Gmsh size limits to 0.20. The second allows Gmsh
to choose characteristic sizes in the 0.10 to 0.20 range. These are meshing
constraints, not exact tetrahedron edge lengths or an exact element count.

Sphere refinement examples
--------------------------

Use the same STL input for each run and change only the size:

.. code-block:: console

   python grid/stl_to_gmsh/generate_volume_mesh.py --stl grid/stl_to_gmsh/3d/sphere/automatic_tetrahedra/python/input/model.stl --output "$HOME/.oneflow/mesh-generation/sphere_coarse.msh" --mesh-size 0.40
   python grid/stl_to_gmsh/generate_volume_mesh.py --stl grid/stl_to_gmsh/3d/sphere/automatic_tetrahedra/python/input/model.stl --output "$HOME/.oneflow/mesh-generation/sphere_medium.msh" --mesh-size 0.20
   python grid/stl_to_gmsh/generate_volume_mesh.py --stl grid/stl_to_gmsh/3d/sphere/automatic_tetrahedra/python/input/model.stl --output "$HOME/.oneflow/mesh-generation/sphere_fine.msh" --mesh-size 0.10

Outputs
-------

For ``--output "$HOME/.oneflow/mesh-generation/sphere.msh"``, the command writes:

* ``$HOME/.oneflow/mesh-generation/sphere.msh`` — MSH 4.1 tetrahedral volume mesh.
* ``$HOME/.oneflow/mesh-generation/sphere.geo`` — generated Gmsh geometry script.
* ``$HOME/.oneflow/mesh-generation/sphere.report.json`` — input topology, mesh counts, Physical
  Groups and validation summary.

The command validates one connected, closed, consistently oriented,
genus-zero STL shell and checks that the generated mesh contains TETRA_4
volume cells, TRI_3 boundary facets, and the requested Physical Groups.
Default group names are ``wall`` (surface) and ``fluid_domain`` (volume);
override them with ``--boundary-name`` and ``--volume-name``.

This command generates a volume inside the closed STL shell. It does not
create an external-flow domain around the STL body. For a full OneFLOW Case
and CFD smoke workflow, see :doc:`stl_to_gmsh_sphere`.
