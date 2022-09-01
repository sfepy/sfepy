Useful Code Snippets and FAQ
============================

.. only:: html

   .. contents:: Table of Contents
      :local:
      :backlinks: top

Code examples below that use `sfepy-*` scripts assume the sfepy package to be
installed, see also :ref:`introduction_installation`.

Miscellaneous
-------------

#. ``No module named 'sfepy.discrete.common.extmods.mappings'``.

   When installing SfePy from sources or using the git version, its extension
   modules have to be compiled before using the package, see
   :ref:`compilation`.

#. The extension modules are compiled in place, but ``ModuleNotFoundError: No
   module named 'sfepy'`` shows up when running some interactive
   examples/scripts/modules from the SfePy source directory.

   On some platforms the current directory is not in the ``sys.path`` directory
   list. Add it using::

     export PYTHONPATH=.

   or add the following code prior to ``sfepy`` imports into the module::

     import sys
     sys.path.append('.')

#. Finite element approximation (field) order and numerical quadrature order.

   SfePy supports reading only straight-facet (linear approximation) meshes,
   nevertheless field orders higher than one can be used, because internally,
   the mesh elements are enriched with the required additional nodes. The
   calculation then occurs on such an augmented mesh with appropriate higher
   order elements.

   The quadrature order equal to two-times the field order (used in many
   examples) works well for bilinear forms with constant (on each element)
   material parameters. For example, a dot product involves integrating ``u *
   v``, so if the approximation order of ``u`` and ``v`` is 1, their product's
   order is 2. Of course, there are terms that could use a lower quadrature
   order, or higher, depending on the data. Increased quadrature order is
   required e.g. in terms with highly oscillating material coefficients.

   Example::

     approx_order = 2
     # The finite element approximation order.
     fields = {
         'displacement': ('real', 3, 'Omega', approx_order),
     }
     # The numerical quadrature order.
     integrals = {
         'i' : 2 * approx_order,
     }


#. Higher order DOF visualization when using an approximation order greater
   than one.

   By default, the additional, higher order DOFs, are not used in the VTK/HDF5
   results files (``'strip'`` linearization kind). To see the influence of
   those DOFs, ``'adaptive'`` linearization has to be used, see
   :ref:`diffusion-sinbc` (declarative API) and
   :ref:`diffusion-laplace_refine_interactive` or
   :ref:`multi_physics-biot_parallel_interactive` (imperative API, search
   ``linearization``).

#. Numbering of DOFs.

   Locally (in a connectivity row), the DOFs are stored DOF-by-DOF (``u_0`` in
   all local nodes, ``u_1`` in all local nodes, ...).

   Globally (in a state vector), the DOFs are stored node-by-node (``u_0, u_1,
   ..., u_X`` in node 0, ``u_0, u_1, ..., u_X`` in node 1, ...).

   See also :func:`create_adof_conn()
   <sfepy.discrete.variables.create_adof_conn>`.

#. Visualization of various FEM-related information.

   - Quadrature rules::

       python3 sfepy/scripts/plot_quadratures.py

   - Facet orientations - run in the source code directory and make sure the
     current directory is in the Python's path list (see
     `Miscellaneous`_)::

       python3 sfepy/postprocess/plot_facets.py

   - Global and local numberings of mesh topological entities (cells, faces,
     edges, vertices)::

       python3 sfepy/scripts/plot_mesh.py meshes/elements/2_4_2.mesh

     The global numbers serve as indices into connectivities. In the plot, the
     global numbers are on the entities, the cell-local ones are inside the
     cells next to each entity towards the cell centroids.

#. How to work with solvers/preconditioners?

   See :ref:`multi_physics-biot_short_syntax` (user-defined preconditioners)
   or :ref:`navier_stokes-stokes_slip_bc` (petsc solver setup).

#. How to get the linear system components: the matrix and the right-hand side?

   To get the residual vector ``r`` (see :doc:`ebcs_implementation`) and the
   tangent matrix ``K``, the imperative API can be used as follows::

     # pb is a Problem instance,
     pb.set_bcs(ebcs=Conditions([...])) # Set Dirichlet boundary conditions.
     pb.set_ics(Conditions([...])) # Set initial conditions (if any).
     variables = pb.get_initial_state()
     pb.time_update()
     pb.update_materials()
     variables.apply_ebc()
     r = pb.equations.eval_residuals(variables())
     K = pb.equations.eval_tangent_matrices(variables(), pb.mtx_a)

   See also :ref:`diffusion-poisson_parallel_interactive`.

#. Where is the code that calculates the element (e.g. stiffness) matrix?

   The code that computes the per element residuals and matrices is organized
   in terms, see :ref:`term_overview` - click on the term class name and then
   "source" link to see the code. The original terms are implemented in C,
   newer terms tend to be implemented directly in Python. The structure and
   attributes of a term class are described in
   :ref:`how_to_implement_a_new_term`.

#. What structural elements (beams, shells, etc.) are available in SfePy?

   The code is currently focused on solid elements. The only supported
   structural element is shell10x, see
   :ref:`linear_elasticity-shell10x_cantilever`.

Mesh-Related Tasks
------------------

#. Checking and fixing a mesh (double vertices, disconnected components, etc.).

   - Show the mesh Euler characteristic, number of components and other
     information::

       sfepy-mesh info -d cylinder.mesh

   - Fix double/disconnected vertices::

       sfepy-convert -m bad.mesh maybe-good.mesh

#. Convert a mesh to another format (as supported by meshio).

   - Simple conversion::

       sfepy-convert mesh.format1 mesh.format2

   - Scaling the mesh anisotropically::

       sfepy-convert -s 2,4,3 cylinder.mesh cylinder-scaled.mesh

#. Verify that regions are correctly defined.

   - Using the problem description files (declarative API)::

       sfepy-run sfepy/examples/diffusion/poisson_short_syntax.py --save-regions-as-groups --solve-not
       sfepy-view -e cylinder_regions.vtk

   - In a script (imperative API)::

       problem.save_regions_as_groups('regions')

#. Remove lower-dimensional entities from a mesh (e.g. edges).

   Use ``sfepy-convert`` with the ``-d <dimension>`` option, where
   ``<dimension>`` is the topological dimension of cells that should be in the
   mesh. For example, ``-d 2`` stores only the 2D cells.

#. It is suggested to use ``msh22`` format instead of the default ``msh4``
   when generating a mesh with ``gmsh``::

      gmsh -2 cylinder.geo -o cylinder.msh -format msh22

   ``msh22`` seems to be more reliable and foolproof when converting.


Regions
-------

#. How to define a region using a function of coordinates in the interactive mode
   (imperative API)?

   Examples:

   -  A facet region defined using a function of mesh vertex coordinates::

       from sfepy.discrete import Function, Functions

       def _get_region(coors, domain=None):
           ii = np.nonzero(coors[:,0] < 0.5)[0]
           return ii

       get_region = Function('get_region', _get_region)
       region = domain.create_region(
           'Region', 'vertices by get_region', 'facet',
           functions=Functions([get_region]),
       )

   - Analogously a cell region defined using the coordinates of cell centroids::

       # ...
       region = domain.create_region(
           'Region', 'cells by get_region', 'cell',
           functions=Functions([get_region]),
       )


Material Parameters
-------------------

#. How to set material parameters per region in the interactive mode
   (imperative API)?

   Example: define ``rho``, ``D`` to have different values in regions ``omega1``,
   ``omega2``::

     m = Material('m', values={'rho': {'omega1': 2700, 'omega2': 6000},
                               'D': {'omega1': D1, 'omega2': D2}})

#. How to implement state dependent materials?

   Besides writing a custom solver, one can use pseudo-time-stepping for this
   purpose, as demonstrated in :ref:`linear_elasticity-material_nonlinearity`
   or :ref:`diffusion-poisson_field_dependent_material`. Note that the examples
   are contrived, and in practice care must be taken to ensure convergence.

#. Why are results of a 2D elasticity simulation not consistent with a properly
   constrained 3D elasticity simulation?

   Possible reason: when using the Young's modulus and Poisson's ratio as input
   parameters, and then calling :func:`stiffness_from_youngpoisson()
   <sfepy.mechanics.matcoefs.stiffness_from_youngpoisson>`, note that the
   default value of the ``plane`` argument is ``'strain'``, corresponding to
   the plane strain assumption, see also :func:`lame_from_youngpoisson()
   <sfepy.mechanics.matcoefs.lame_from_youngpoisson>`. Try setting
   ``plane='stress'``.

#. How to set (time-dependent) material parameters by a function in the
   interactive mode (imperative API)?

   Example (also showing the full material function signature)::

     from sfepy.discrete import Material, Function

     def get_pars(ts, coors, mode=None,
                  equations=None, term=None, problem=None, **kwargs):
         value1 = a_function(ts.t, coors)
         value2 = another_function(ts.step, coors)
         if mode == 'qp':
             out = {
                 'value1' : value1.reshape(coors.shape[0], 1, 1),
                 'value2' : value2.reshape(coors.shape[0], 1, 1),
             }
             return out
     m = Material('m', function=Function('get_pars', get_pars))

#. How to get cells corresponding to coordinates in a material function?

   The full signature of the material function is::

     def get_pars(ts, coors, mode=None,
                  equations=None, term=None, problem=None, **kwargs)

   Thus it has access to ``term.region.cells``, hence access to the cells that
   correspond to the coordinates. The length of the ``coors`` is
   ``n_cell * n_qp``, where ``n_qp`` is the number of quadrature
   points per cell, and ``n_cell = len(term.region.cells)``, so that
   ``coors.reshape((n_cell, n_qp, -1))`` can be used.
