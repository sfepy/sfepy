Tutorial
========

*SfePy* can be used in two basic ways:
  #. a black-box partial differential equation (PDE) solver,
  #. a Python package to build custom applications involving solving PDEs by the
     finite element (FE) method.

This tutorial focuses on the first way and introduces the basic concepts and
nomenclature used in the following parts of the documentation.

Notes on solving PDEs by the Finite Element Method
--------------------------------------------------

Topics which should eventually be discussed here
  * A description of the PDE weak form
  * Discussion of discretization and meshing

It is planned to have an example based on the Poisson's equation here. For now,
please refer to the wikipedia page at
http://en.wikipedia.org/wiki/Finite_element_method for a basic description of
the Finite Element Method.

Running a simulation
--------------------

The following commands should be run in the top-level directory of the *SfePy*
source tree after compiling the C extension files. See
:ref:`introduction_installation` for full installation instructions. 

Running *SfePy* through the GUI
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The easisest way to run *SfePy* is through the GUI utility. 

* Locate the ``sfepy_gui.py`` script in the top-level directory

* Most environments allow you to double-click on this file to execute the
  python script. Consult your platform documentation to find out how to
  associate python scripts with the python interpreter.

* After executing the ``sfepy_gui.py`` script, you should see the following
  window:

.. image:: images/sfepy_gui.png

* Click the *Browse* button to the right of the *input file name* box and
  browse for ``examples/diffusion/poisson.py`` as indicated in the figure

* The solution will start automatically and an output box will pop up when the
  simulation completes

* Viewing the output is described in the :ref:`postprocessing` section below

Invoking *SfePy* from the command line
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This section introduces the basics of running *SfePy* on the command line. The
``$`` indicates the command prompt of your terminal.

* The script ``simple.py`` is the most basic starting point in *SfePy*. It is
  invoked as follows::

    $ ./simple.py examples/diffusion/poisson.py

  * ``examples/diffusion/poisson.py`` is the *SfePy* *problem description*
    file, which defines the problem to be solved in terms *SfePy* can
    understand

  * Running the above command creates the output file ``cylinder.vtk`` in the
    *SfePy* top-level directory

* *SfePy* can also be invoked interactively with the ``isfepy`` script::

    $ ./isfepy

  * Follow the help information printed on startup to solve the
    Poisson's equation example above

.. _postprocessing:

Postprocessing the results
^^^^^^^^^^^^^^^^^^^^^^^^^^

* The ``postproc.py`` script can be used for quick postprocessing and
  visualization of the *SfePy* output files. It requires mayavi2 installed on
  your system.

  * As a simple example, try::

    $ ./postproc.py cylinder.vtk

  * The following interactive 3D window should display:

.. image:: images/postproc_simple.png

* The left mouse button by itself orbits the 3D view

* Holding shift and the left mouse button pans the view

* Holding control and the left mouse button rotates about the screen normal axis

* The right mouse button controls the zoom


Example problem description file
--------------------------------

Here we discuss the contents of the ``examples/diffusion/poisson.py`` problem
description file. For additional examples, see the problem description files in
the ``examples/`` directory of SfePy.

Open the ``examples/diffusion/poisson.py`` file in your favorite text
editor. Note that the file is a regular python source code.

:: 

    #! Poisson Equation
    #! ================
    #$ \centerline{Example input file, \today}

    #! Mesh
    #! ----
    from sfepy import data_dir

    filename_mesh = data_dir + '/meshes/3d/cylinder.mesh'

The ``filename_mesh`` variable points to the file containing the mesh for the
particular problem. *SfePy* supports a variety of mesh formats.

::

    #! Materials
    #! ---------
    #$ Here we define just a constant coefficient $c$ of the Poisson equation,
    #$ using the 'values' attribute. Other possible attribute is 'function', for
    #$ material coefficients computed/obtained at runtime.

    material_2 = {
        'name' : 'coef',
        'values' : {'val' : 1.0},
    }

Many finite element problems require the definition of material parameters.
These can be handled in *SfePy* with material variables which associate the
material parameters with the corresponding region of the mesh.

::

    #! Fields
    #! ------
    #! A field is used mainly to define the approximation on a (sub)domain, i.e. to
    #$ define the discrete spaces $V_h$, where we seek the solution.
    #!
    #! The Poisson equation can be used to compute e.g. a temperature distribution,
    #! so let us call our field 'temperature'. On a region called 'Omega'
    #! (see below) it will be approximated using P1 finite elements.

    field_1 = {
        'name' : 'temperature',
        'dim' : (1,1),
        'domain' : 'Omega',
        'bases' : {'Omega' : '3_4_P1'}
    }



::

    #! Variables
    #! ---------
    #! One field can be used to generate discrete degrees of freedom (DOFs) of
    #! several variables. Here the unknown variable (the temperature) is called
    #! 't', it's associated DOF name is 't.0' --- this will be referred to
    #! in the Dirichlet boundary section (ebc). The corresponding test variable of
    #! the weak formulation is called 's'. Notice that the 'dual' item of a test
    #! variable must specify the unknown it corresponds to.

    variable_1 = {
        'name' : 't',
        'kind' : 'unknown field',
        'field' : 'temperature',
        'order' : 0, # order in the global vector of unknowns
    }

    variable_2 = {
        'name' : 's',
        'kind' : 'test field',
        'field' : 'temperature',
        'dual' : 't',
    }



::

    #! Regions
    #! -------
    region_1000 = {
        'name' : 'Omega',
        'select' : 'elements of group 6',
    }

    region_03 = {
        'name' : 'Gamma_Left',
        'select' : 'nodes in (x < 0.00001)',
    }

    region_4 = {
        'name' : 'Gamma_Right',
        'select' : 'nodes in (x > 0.099999)',
    }

Regions assign names to various parts of the finite element mesh. The region
names can later be referred to, for example when specifying portions of the mesh
to apply boundary conditions to. Regions can be specified in a variety of ways,
including by element or by node. Here, Omega is the elemental domain over which
the PDE is solved and Gamma_Left and Gamma_Right define surfaces upon which the
boundary conditions will be applied.

::

    #! Boundary Conditions
    #! -------------------
    #! Essential (Dirichlet) boundary conditions can be specified as follows:
    ebc_1 = {
        'name' : 't1',
        'region' : 'Gamma_Left',
        'dofs' : {'t.0' : 2.0},
    }

    ebc_2 = {
        'name' : 't2',
        'region' : 'Gamma_Right',
        'dofs' : {'t.0' : -2.0},
    }

Boundary conditions place restrictions on the finite element formulation and
create a unique solution to the problem. Here, we specify that a temperature of
+2 is applied to the left surface of the mesh and a temperature of -2 is applied
to the right surface.

::

    #! Equations
    #! ---------
    #$ The weak formulation of the Poisson equation is:
    #$ \begin{center}
    #$ Find $t \in V$, such that
    #$ $\int_{\Omega} c\ \nabla t : \nabla s = f, \quad \forall s \in V_0$.
    #$ \end{center}
    #$ The equation below directly corresponds to the discrete version of the
    #$ above, namely:
    #$ \begin{center}
    #$ Find $\bm{t} \in V_h$, such that
    #$ $\bm{s}^T (\int_{\Omega_h} c\ \bm{G}^T G) \bm{t} = 0, \quad \forall \bm{s}
    #$ \in V_{h0}$,
    #$ \end{center}
    #$ where $\nabla u \approx \bm{G} \bm{u}$. Below we use $f = 0$ (Laplace
    #$ equation).
    #! We also define an integral here: 'gauss_o1_d3' says that we wish to use
    #! quadrature of the first order in three space dimensions.
    integral_1 = {
        'name' : 'i1',
        'kind' : 'v',
        'quadrature' : 'gauss_o2_d3',
    }

Integrals specify which numerical scheme to use. Here we are using a 2nd order
quadrature over a 3 dimensional space.

::

    equations = {
        'Temperature' : """dw_laplace.i1.Omega( coef.val, s, t ) = 0"""
    }

The equations block is the heart of the *SfePy* problem definition file. Here,
we are specifying that the Laplacian of the temperature (in the weak
formulation) is 0, where ``coef.val`` is a material constant. We are using the
``i1`` integral defined previously, over the domain specified by the region
Omega.

The above syntax is useful for defining *custom integrals* with
user-defined quadrature points and weights, see :ref:`ug_integrals`. The
above uniform integration can be more easily achieved by::

    equations = {
        'Temperature' : """dw_laplace.2.Omega( coef.val, s, t ) = 0"""
    }

The integration order is specified directly in place of the integral
name. The integral definition is superfluous in this case.

::

    #! Linear solver parameters
    #! ---------------------------
    #! Use umfpack, if available, otherwise superlu.
    solver_0 = {
        'name' : 'ls',
        'kind' : 'ls.scipy_direct',
        'method' : 'auto',
    }

Here, we specify which kind of solver to use for linear equations.

::

    #! Nonlinear solver parameters
    #! ---------------------------
    #! Even linear problems are solved by a nonlinear solver (KISS rule) - only one
    #! iteration is needed and the final rezidual is obtained for free.
    solver_1 = {
        'name' : 'newton',
        'kind' : 'nls.newton',

        'i_max'      : 1,
        'eps_a'      : 1e-10,
        'eps_r'      : 1.0,
        'macheps'   : 1e-16,
        'lin_red'    : 1e-2, # Linear system error < (eps_a * lin_red).
        'ls_red'     : 0.1,
        'ls_red_warp' : 0.001,
        'ls_on'      : 1.1,
        'ls_min'     : 1e-5,
        'check'     : 0,
        'delta'     : 1e-6,
        'is_plot'    : False,
        'problem'   : 'nonlinear', # 'nonlinear' or 'linear' (ignore i_max)
    }

Here, we specify the nonlinear solver kind and options. The convergence
parameters can be adjusted if necessary, otherwise leave the default.

::

    #! Options
    #! -------
    #! Use them for anything you like... Here we show how to tell which solvers
    #! should be used - reference solvers by their names.
    options = {
        'nls' : 'newton',
        'ls' : 'ls',
    }

The solvers to use are specified in the options block. We can define multiple
solvers with different convergence parameters if necessary.

::

    #! FE assembling parameters
    #! ------------------------
    #! 'chunk_size' is now unused, deprecated, and will be removed.
    fe = {
        'chunk_size' : 1000
    }

Interactive Example: Linear Elasticity
--------------------------------------

This example shows how to use *SfePy* interactively, but also how to
make a custom simulation script. We will use ``isfepy`` for the
explanation, but regular Python shell, or IPython would do as well,
provided the proper modules are imported (see the help information
printed on startup of ``isfepy``).

We wish to solve the following linear elasticity problem:

.. math::
   :label: eq_linear_elasticity

    - \pdiff{\sigma_{ij}(\ul{u})}{x_j} + f_i = 0 \mbox{ in }\Omega,
    \quad \ul{u} = 0 \mbox{ on } \Gamma_1,
    \quad u_1 = \bar{u}_1 \mbox{ on } \Gamma_2 \;,

where the stress is defined as :math:`\sigma_{ij} = 2 \mu e_{ij} +
\lambda e_{kk} \delta_{ij}`, :math:`\lambda`, :math:`\mu` are the Lamé's
constants, the strain is :math:`e_{ij}(\ul{u}) =
\frac{1}{2}(\pdiff{u_i}{x_j} + \pdiff{u_j}{x_i})` and :math:`\ul{f}` are
volume forces. This can be written in general form as
:math:`\sigma_{ij}(\ul{u}) = D_{ijkl} e_{kl}(\ul{u})`, where in our case
:math:`D_{ijkl} = \mu (\delta_{ik} \delta_{jl}+\delta_{il}
\delta_{jk}) + \lambda \ \delta_{ij} \delta_{kl}`.

In the weak form the equation :eq:`eq_linear_elasticity` is

.. math::
   :label: eq_wlinear_elasticity

    \int_{\Omega} D_{ijkl} e_{kl}(\ul{u}) e_{kl}(\ul{v}) + \int_{\Omega}
    f_i v_i = 0 \;,

where :math:`\ul{v}` is the test function, and both :math:`\ul{u}`,
:math:`\ul{v}` belong to a suitable function space.

**Hint:** Whenever you create a new object (e.g. a Mesh instance, see
below), try to print it using the `print` statement - it will give you
insight about the object internals.

The whole example summarized in a script is below in
:ref:`tutorial_interactive_source`.

Run the ``isfepy`` script::

    $ ./isfepy

The output should look like this:

.. sourcecode:: ipython

    Python 2.6.5 console for SfePy 2010.2-git-11cfd34 (8c4664610ed4b85851966326aaa7ce36e560ce7a)

    These commands were executed:
    >>> from sfepy.base.base import *
    >>> from sfepy.fem import *
    >>> from sfepy.applications import pde_solve
    >>> from sfepy.postprocess import Viewer

    When in SfePy source directory, try:
    >>> pb, vec, data = pde_solve('examples/diffusion/poisson.py')
    >>> view = Viewer(pb.get_output_name())
    >>> view()

    When in another directory (and SfePy is installed), try:
    >>> from sfepy import data_dir
    >>> pb, vec, data = pde_solve(data_dir + '/examples/diffusion/poisson.py')
    >>> view = Viewer(pb.get_output_name())
    >>> view()

    Documentation can be found at http://sfepy.org

    In [1]:

Read a finite element mesh, that defines the domain :math:`\Omega`.

.. sourcecode:: ipython

    In [1]: mesh = Mesh.from_file('meshes/2d/rectangle_tri.mesh')

Create a domain. The domain allows defining regions or subdomains.

.. sourcecode:: ipython

    In [2]: domain = Domain('domain', mesh)
    sfepy: setting up domain edges...
    sfepy: ...done in 0.01 s

Define the regions - the whole domain :math:`\Omega`, where the solution
is sought, and :math:`\Gamma_1`, :math:`\Gamma_2`, where the boundary
conditions will be applied. As the domain is rectangular, we first get a
bounding box to get correct bounds for selecting the boundary edges.

.. sourcecode:: ipython

    In [3]: min_x, max_x = domain.get_mesh_bounding_box()[:,0]
    In [4]: eps = 1e-8 * (max_x - min_x)
    In [5]: omega = domain.create_region('Omega', 'all')
    In [6]: gamma1 = domain.create_region('Gamma1',
       ...:                               'nodes in x < %.10f' % (min_x + eps))
    In [7]: gamma2 = domain.create_region('Gamma2',
       ...:                               'nodes in x > %.10f' % (max_x - eps))

Next we define the actual finite element approximation using the
:class:`Field` class.

.. sourcecode:: ipython

    In [8]: field = Field('fu', nm.float64, 'vector', omega,
       ...:               space='H1', poly_space_base='lagrange', approx_order=2)

Using the field `fu`, we can define both the unknown variable :math:`\ub` and
the test variable :math:`\vb`.

.. sourcecode:: ipython

    In [9]: u = FieldVariable('u', 'unknown', field, mesh.dim)
    In [10]: v = FieldVariable('v', 'test', field, mesh.dim,
       ....:                   primary_var_name='u')

Before we can define the terms to build the equation of linear
elasticity, we have to create also the materials, i.e. define the
(constitutive) parameters. The linear elastic material `m` will be
defined using the two Lamé constants :math:`\lambda = 1`, :math:`\mu =
1`. The volume forces will be defined also as a material, as a constant
(column) vector :math:`[0.02, 0.01]^T`.

.. sourcecode:: ipython

    In [11]: m = Material('m', lam=1.0, mu=1.0)
    In [12]: f = Material('f', val=[[0.02], [0.01]])

One more thing needs to be defined - the numerical quadrature that will
be used to integrate each term over its domain.

.. sourcecode:: ipython

    In [14]: integral = Integral('i', order=3)

Now we are ready to define the two terms and build the equations.

.. sourcecode:: ipython

    In [15]: from sfepy.terms import Term
    In [16]: t1 = Term.new('dw_lin_elastic_iso(m.lam, m.mu, v, u)',
                  integral, omega, m=m, v=v, u=u)
    In [17]: t2 = Term.new('dw_volume_lvf(f.val, v)', integral, omega, f=f, v=v)
    In [18]: eq = Equation('balance', t1 + t2)
    In [19]: eqs = Equations([eq])
    sfepy: setting up dof connectivities...
    sfepy: ...done in 0.00 s
    sfepy: describing geometries...
    sfepy: ...done in 0.00 s

The equations have to be completed by boundary conditions. Let us clamp
the left edge :math:`\Gamma_1`, and shift the right edge
:math:`\Gamma_2` in the :math:`x` direction a bit, depending on the
:math:`y` coordinate.

.. sourcecode:: ipython

   In [20]: from sfepy.fem.conditions import Conditions, EssentialBC
   In [21]: fix_u = EssentialBC('fix_u', gamma1, {'u.all' : 0.0})
   In [22]: def shift_u_fun(ts, coors, bc=None, shift=0.0):
      ....:     val = shift * coors[:,1]**2
      ....:     return val
      ....:
   In [23]: bc_fun = Function('shift_u_fun', shift_u_fun,
      ....:                   extra_args={'shift' : 0.01})
   In [24]: shift_u = EssentialBC('shift_u', gamma2, {'u.0' : bc_fun})

The last thing to define before building the problem are the
solvers. Here we just use a sparse direct SciPy solver and the SfePy
Newton solver with default parameters. We also wish to store the
convergence statistics of the Newton solver. As the problem is linear,
it should converge in one iteration.

.. sourcecode:: ipython

    In [25]: from sfepy.solvers.ls import ScipyDirect
    In [26]: from sfepy.solvers.nls import Newton
    In [27]: ls = ScipyDirect({})
    In [28]: nls_status = IndexedStruct()
    In [29]: nls = Newton({}, lin_solver=ls, status=nls_status)

Now we are ready to create a :class:`ProblemDefinition` instance. Note
that the step above is not really necessary - the above solvers are
constructed by default. We did them to get the `nls_status`.

.. sourcecode:: ipython

    In [30]: pb = ProblemDefinition('elasticity', equations=eqs, nls=nls, ls=ls)

The :class:`ProblemDefinition` has several handy methods for
debugging. Let us try saving the regions into a VTK file.

.. sourcecode:: ipython

    In [31]: pb.save_regions_as_groups('regions')
    sfepy: saving regions as groups...
    sfepy:   Omega
    sfepy:   Gamma1
    sfepy:   Gamma2
    sfepy:   Gamma1
    sfepy: ...done

And view them.

.. sourcecode:: ipython

    In [32]: view = Viewer('regions.vtk')
    In [33]: view()
    sfepy: point scalars Gamma1 [ 0.  0.  0.]
    sfepy: point scalars Gamma2 [ 11.   0.   0.]
    sfepy: point scalars Omega [ 22.   0.   0.]
    Out[33]: <sfepy.postprocess.viewer.ViewerGUI object at 0x93ea5f0>

You should see this:

.. image:: images/linear_elasticity_regions.png

Finally, we apply the boundary conditions, solve the problem, save and
view the results.

.. sourcecode:: ipython

    In [34]: pb.time_update(ebcs=Conditions([fix_u, shift_u]))
    sfepy: updating materials...
    sfepy:     m
    sfepy:     f
    sfepy: ...done in 0.01 s
    sfepy: updating variables...
    sfepy: ...done
    sfepy: matrix shape: (1815, 1815)
    sfepy: assembling matrix graph...
    sfepy: ...done in 0.00 s
    sfepy: matrix structural nonzeros: 39145 (1.19e-02% fill)
    In [35]: vec = pb.solve()
    sfepy: nls: iter: 0, residual: 1.343114e+01 (rel: 1.000000e+00)
    sfepy:   rezidual:    0.00 [s]
    sfepy:      solve:    0.01 [s]
    sfepy:     matrix:    0.00 [s]
    sfepy: nls: iter: 1, residual: 2.567997e-14 (rel: 1.911972e-15)
    In [36]: print nls_status
    -------> print(nls_status)
    IndexedStruct
      condition:
        0
      err:
        2.56799662867e-14
      err0:
        13.4311385972
      time_stats:
        {'rezidual': 0.0, 'solve': 0.010000000000001563, 'matrix': 0.0}
    In [37]: pb.save_state('linear_elasticity.vtk', vec)
    In [38]: view = Viewer('linear_elasticity.vtk')
    In [39]: view()
    sfepy: point vectors u [ 0.  0.  0.]
    Out[39]: <sfepy.postprocess.viewer.ViewerGUI object at 0xad61bf0>

This is the resulting image:

.. image:: images/linear_elasticity_solution1.png

The default view is not very fancy. Let us show the displacements by
shifting the mesh. Close the previous window and do:

.. sourcecode:: ipython

    In [56]: view(vector_mode='warp_norm', rel_scaling=2,
       ....:      is_scalar_bar=True, is_wireframe=True)
    sfepy: point vectors u [ 0.  0.  0.]
    Out[56]: <sfepy.postprocess.viewer.ViewerGUI object at 0xad61bf0>

And the result is:

.. image:: images/linear_elasticity_solution2.png

See the docstring of `view()` and play with its options.

.. _tutorial_interactive_source:

Complete Example as a Script
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The source code: :download:`linear_elasticity.py
<../examples/standalone/interactive/linear_elasticity.py>`. It should be
run from the *SfePy* source directory so that it finds the mesh file.

.. literalinclude:: ../examples/standalone/interactive/linear_elasticity.py
   :linenos:

