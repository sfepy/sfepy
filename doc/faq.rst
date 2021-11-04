Useful Code Snippets and F.A.Q.
-------------------------------

#. How to set material parameters per region in the interactive mode
   (imperative API)?

Example: define `rho`, `D` to have different values in regions `omega1`,
`omega2`::

  m = Material('m', values={'rho': {'omega1': 2700, 'omega2': 6000},
                            'D': {'omega1': D1, 'omega2': D2}})

Mesh-related Tasks
^^^^^^^^^^^^^^^^^^

#. Checking and fixing a mesh (double vertices, disconnected components, etc.).

- Show the mesh Euler characteristic, number of components and other
  information::

    python3 script/show_mesh_info.py -d cylinder.mesh

- Fix double/disconnected vertices::

    python3 script/convert_mesh.py -m bad.mesh maybe-good.mesh

#. Convert a mesh to another format (as supported by meshio).

- Simple conversion::

    python3 script/convert_mesh.py mesh.format1 mesh.format2

- Scaling the mesh anisotropically::

    python3 script/convert_mesh.py -s 2,4,3 cylinder.mesh cylinder-scaled.mesh
