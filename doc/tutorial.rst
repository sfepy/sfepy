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
    #! 'chunk_size' determines maximum number of elements to assemble in one C
    #! function call. Higher values mean faster assembling, but also more memory
    #! usage.
    fe = {
        'chunk_size' : 1000
    }

The ``chunk_size`` parameter can be used to tweak the tradeoff between faster
CPU and higher memory usage.
