.. include:: links.inc

.. _sec-tutorial:

Tutorial
========

.. contents:: Table of Contents
   :local:
   :backlinks: top

*SfePy* can be used in two basic ways:
  #. a black-box partial differential equation (PDE) solver,
  #. a Python package to build custom applications involving solving PDEs by the
     finite element (FE) method.

This tutorial focuses on the first way and introduces the basic concepts
and nomenclature used in the following parts of the documentation. Check
also the :doc:`primer` which focuses on a particular problem in detail.

Notes on solving PDEs by the Finite Element Method
--------------------------------------------------

The Finite Element Method (FEM) is the numerical method for solving Partial
Differential Equations (PDEs). FEM was developed in the middle of XX. century and
now it is widely used in different areas of science and engineering, including
mechanical and structural design, biomedicine, electrical and power design,
fluid dynamics and other. FEM is based on a very elegant mathematical theory of
weak solution of PDEs. In this section we will briefly discuss basic ideas
underlying FEM.

.. _poisson-weak-form-tutorial:

Strong form of Poisson's equation and its integration
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Let us start our discussion about FEM with the strong form of Poisson's
equation

.. math::
   :label: eq_spoisson

   \Delta T = f(x), \quad x \in \Omega,

.. math::
   :label: eq_spoisson_dbc

   T = u(x), \quad x \in \Gamma_D,

.. math::
   :label: eq_spoisson_nbc

   \nabla T \cdot \mathbf{n} = g(x), \quad x \in \Gamma_N,

where :math:`\Omega \subset \mathbb{R}^n` is the solution domain with the
boundary :math:`\partial \Omega`, :math:`\Gamma_D` is the part of the boundary
where Dirichlet boundary conditions are given, :math:`\Gamma_N` is the part of
the boundary where Neumann boundary conditions are given, :math:`T(x)` is the
unknown function to be found, :math:`f(x), u(x), g(x)` are known functions.

FEM is based on a weak formulation. The weak form of the
equation :eq:`eq_spoisson` is

.. math::
    \int\limits_{\Omega} (\Delta T - f) \cdot s \, \mathrm{d}\Omega = 0,

where :math:`s` is a **test** function. Integrating this equation by parts

.. math::
  0 = \int\limits_{\Omega} (\Delta T - f) \cdot s \, \mathrm{d}\Omega
    = \int\limits_{\Omega} \nabla \cdot (\nabla T) \cdot s \, \mathrm{d}\Omega
      - \int_{\Omega} f \cdot s \, \mathrm{d}\Omega =

.. math::
    = - \int\limits_{\Omega} \nabla T \cdot \nabla s \, \mathrm{d}\Omega
      + \int\limits_{\Omega} \nabla \cdot (\nabla T \cdot s) \, \mathrm{d}\Omega
      - \int\limits_{\Omega} f \cdot s \, \mathrm{d}\Omega

and applying Gauss theorem we obtain:

.. math::
  0 = - \int\limits_{\Omega} \nabla T \cdot \nabla s \, \mathrm{d}\Omega
      + \int\limits_{\Gamma_D \cup \Gamma_N} \!\!\!\!
            s \cdot (\nabla T \cdot \mathbf{n}) \, \mathrm{d}\Gamma
      - \int\limits_{\Omega} f \cdot s \, \mathrm{d}\Omega

or

.. math::
    \int\limits_{\Omega} \nabla T \cdot \nabla s \, \mathrm{d}\Omega
  =
    \int\limits_{\Gamma_D \cup \Gamma_N} \!\!\!\!
        s \cdot (\nabla T \cdot \mathbf{n}) \, \mathrm{d}\Gamma
  - \int\limits_{\Omega} f \cdot s \, \mathrm{d}\Omega.

The surface integral term can be split into two integrals, one over the
Dirichlet part of the surface and second over the Neumann part

.. math::
    :label: eq_wpoisson_full

    \int\limits_{\Omega} \nabla T \cdot \nabla s \, \mathrm{d}\Omega
    =
    \int\limits_{\Gamma_D}
        s \cdot (\nabla T \cdot \mathbf{n}) \, \mathrm{d}\Gamma
    + \int\limits_{\Gamma_N}
        s \cdot (\nabla T \cdot \mathbf{n}) \, \mathrm{d}\Gamma
    - \int\limits_{\Omega} f \cdot s \, \mathrm{d}\Omega.

The equation :eq:`eq_wpoisson_full` is the initial weak form of the
Poisson's problem :eq:`eq_spoisson`--:eq:`eq_spoisson_nbc`. But we can not work
with it without applying the boundary conditions. So it is time to talk about
the boundary conditions.

Dirichlet Boundary Conditions
"""""""""""""""""""""""""""""

On the Dirichlet part of the surface we have two restrictions. One is the Dirichlet
boundary conditions :math:`T(x) = u(x)` as they are, and the second is the
integral term over :math:`\Gamma_D` in equation :eq:`eq_wpoisson_full`. To be
consistent we have to use only the Dirichlet conditions and avoid the integral
term. To implement this we can take the function :math:`T \in V(\Omega)` and
the test function :math:`s \in V_0(\Omega)`, where

.. math::
    V(\Omega) = \{f(x) \in H^1(\Omega)\},

.. math::
  V_0(\Omega) = \{f(x) \in H^1(\Omega); f(x) = 0, x \in \Gamma_D\}.

In other words the unknown function :math:`T` must be continuous together with
its gradient in the domain. In contrast the test function :math:`s` must be
also continuous together with its gradient in the domain but it should be zero
on the surface :math:`\Gamma_D`.

With this requirement the integral term over Dirichlet part of the surface
is vanishing and the weak form of the Poisson equation for
:math:`T \in V(\Omega)` and :math:`s \in V_0(\Omega)` becomes

.. math::
    \int\limits_{\Omega} \nabla T \cdot \nabla s \, \mathrm{d}\Omega
    =
    \int\limits_{\Gamma_N}
        s \cdot (\nabla T \cdot \mathbf{n}) \, \mathrm{d}\Gamma
    - \int\limits_{\Omega} f \cdot s \, \mathrm{d}\Omega,

    T(x) = u(x), \quad x \in \Gamma_D.

That is why Dirichlet conditions in FEM terminology are called
**Essential Boundary Conditions**. These conditions are not a part of the weak
form and they are used as they are.

Neumann Boundary Conditions
"""""""""""""""""""""""""""

The Neumann boundary conditions correspond to the known flux
:math:`g(x) = \nabla T \cdot \mathbf{n}`. The integral term over the Neumann
surface in the equation :eq:`eq_wpoisson_full` contains exactly the same flux.
So we can use the known function :math:`g(x)` in the integral term:

.. math::
    \int\limits_{\Omega} \nabla T \cdot \nabla s \, \mathrm{d}\Omega
    =
    \int\limits_{\Gamma_N} g \cdot s \, \mathrm{d}\Gamma
    - \int\limits_{\Omega} f \cdot s \, \mathrm{d}\Omega,

where test function :math:`s` also belongs to the space :math:`V_0`.

That is why Neumann conditions in FEM terminology are called
**Natural Boundary Conditions**. These conditions are a part of weak form terms.

The weak form of the Poisson's equation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Now we can write the resulting weak form for the Poisson's problem
:eq:`eq_spoisson`--:eq:`eq_spoisson_nbc`. For any test function
:math:`s \in V_0(\Omega)` find :math:`T \in V(\Omega)`
such that

.. math::
    :label: eq_wpoisson_final

    \boxed {
    \begin{split}
      \int\limits_{\Omega} \nabla T \cdot \nabla s \, \mathrm{d}\Omega
      & =
      \int\limits_{\Gamma_N} g \cdot s \, \mathrm{d}\Gamma
      - \int\limits_{\Omega} f \cdot s \, \mathrm{d}\Omega, \quad \mbox{and}\\
      T(x) & = u(x), \quad x \in \Gamma_D.
    \end{split}
    }

Discussion of discretization and meshing
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

It is planned to have an example of the discretization based on the Poisson's
equation weak form :eq:`eq_wpoisson_final`. For now, please refer to the
wikipedia page `Finite Element Method`_ for a basic description of the
disretization and meshing.


Numerical solution of the problem
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To solve numerically given problem based on the weak form
:eq:`eq_wpoisson_final` we have to go through 5 steps:

#. Define geometry of the domain :math:`\Omega` and
   surfaces :math:`\Gamma_D` and :math:`\Gamma_N`.
#. Define the known functions :math:`f`, :math:`u` and :math:`g`.
#. Define the unknown function :math:`T` and the test functions :math:`s`.
#. Define essential boundary conditions (Dirichlet conditions)
   :math:`T(x) = u(x), x \in \Gamma_D.`
#. Define equation and natural boundary conditions (Neumann conditions)
   as the set of all integral terms
   :math:`\int\limits_{\Omega} \nabla T \cdot \nabla s \, \mathrm{d}\Omega`,
   :math:`\int\limits_{\Gamma_N} g \cdot s \, \mathrm{d}\Gamma`,
   :math:`\int\limits_{\Omega} f \cdot s \, \mathrm{d}\Omega`.

In the next section we will discuss how to define all these things in SfePy.

Basic notions
-------------

The simplest way of using *SfePy* is to solve a system of PDEs defined
in a **problem description file**, also referred to as **input
file**. In such a file, the problem is described using several keywords
that allow one to define the equations, variables, finite element
approximations, solvers, solution domain and subdomains etc., see
:ref:`sec-problem-description-file` for a full list of those keywords.

The syntax of the problem description file is very simple yet powerful,
as the file itself is just a regular Python module that can be normally
imported - no special parsing is necessary. The keywords mentioned above
are regular Python variables (usually of the `dict` type) with special
names. Historically, the keywords exist in two flavors:

* **long syntax** is the original one - it is longer to type, but the
  individual fields are named, so it might be easier/understand to read
  for newcomers.
* **short syntax** was added later to offer brevity for "expert" use.

Below we show:

#. how to solve a problem given by a problem description file, and
#. explain the elements of the file on several examples.

But let us begin with a slight detour...

Sneak peek: what is going on under the hood
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

#. A top-level script (usually ``simple.py``, as in this tutorial) reads
   in an input file.

#. Following the contents of the input file, a :class:`ProblemDefinition
   <sfepy.fem.problemDef.ProblemDefinition>` instance is created - this
   is the input file coming to life. Let us call the instance
   ``problem``.

   * The ``problem`` sets up its domain, regions (various sub-domains),
     fields (the FE approximations), the equations and the solvers. The
     equations determine the materials and variables in use - only those
     are fully instantiated, so the input file can safely contain
     definitions of items that are not used actually.

#. Prior to solution, ``problem.time_update()`` function has to be
   called to setup boundary conditions, material parameters and other
   potentially time-dependent data. This holds also for stationary
   problems with a single "time step".

#. The solution is then obtained by calling ``problem.solve()``
   function.

#. Finally, the solution can be stored using ``problem.save_state()``

The above last three steps are essentially repeated for each time
step. So that is it - using the code a black-box PDE solver shields the
user from having to create the :class:`ProblemDefinition
<sfepy.fem.problemDef.ProblemDefinition>` instance by hand. But note
that this is possible, and often necessary when the flexibility of the
default solvers is not enough. At the end of the tutorial an example
demonstrating the interactive creation of the ``problem`` is shown, see
:ref:`sec-interactive-example-linear-elasticity`.

Now let us continue with running a simulation.

Running a simulation
--------------------

The following commands should be run in the top-level directory of the *SfePy*
source tree after compiling the C extension files. See
:ref:`introduction_installation` for full installation instructions.

.. _invoking_command_line:

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
   :width: 70 %
   :align: center

* The left mouse button by itself orbits the 3D view

* Holding shift and the left mouse button pans the view

* Holding control and the left mouse button rotates about the screen normal axis

* The right mouse button controls the zoom

.. _poisson-example-tutorial:

Example problem description file
--------------------------------

Here we discuss the contents of the :download:`examples/diffusion/poisson.py
<../examples/diffusion/poisson.py>` problem
description file. For additional examples, see the problem description files in
the ``examples/`` directory of *SfePy*.

The problem at hand is the following:

.. math::
   :label: eq_laplace_static

    c \Delta T = f \mbox{ in }\Omega,\quad T(t) = \bar{T}(t)
    \mbox{ on } \Gamma \;,

where :math:`\Gamma \subseteq \Omega` is a subset of the domain :math:`\Omega`
boundary. For simplicity, we set :math:`f \equiv 0`, but we still work with
the material constant :math:`c` even though it has no influence on the
solution in this case. We also assume zero fluxes over :math:`\partial
\Omega \setminus \Gamma`, i.e. :math:`\pdiff{T}{\ul{n}} = 0` there.  The
particular boundary conditions used below are :math:`T = 2` on the left side
of the cylindrical domain depicted in the previous section and :math:`T = -2`
on the right side.

The first step to do is to write :eq:`eq_laplace_static` in *weak
formulation* :eq:`eq_wpoisson_final`. The :math:`f = 0`,
:math:`g = \pdiff{T}{\ul{n}} = 0`. So only one term in weak form
:eq:`eq_wpoisson_final` remains:

.. math::
   :label: eq_wlaplace_static

    \int_{\Omega} c\ \nabla T \cdot \nabla s = 0, \quad \forall s \in V_0 \;.

Comparing the above integral term with the
long table in :ref:`term_overview`, we can see that *SfePy* contains
this term under name `dw_laplace`. We are now ready to proceed to the
actual problem definition.

Long syntax of keywords
^^^^^^^^^^^^^^^^^^^^^^^

The example uses **long syntax** of the keywords. In next subsection, we
show the same example written in **short syntax**.

Open the :download:`examples/diffusion/poisson.py
<../examples/diffusion/poisson.py>` file in your favorite text
editor. Note that the file is a regular python source code.

::

    from sfepy import data_dir

    filename_mesh = data_dir + '/meshes/3d/cylinder.mesh'

The ``filename_mesh`` variable points to the file containing the mesh for the
particular problem. *SfePy* supports a variety of mesh formats.

::

    material_2 = {
        'name' : 'coef',
        'values' : {'val' : 1.0},
    }

Here we define just a constant coefficient :math:`c` of the Poisson
equation, using the ``'values'`` attribute. Other possible attribute is
``'function'``, for material coefficients computed/obtained at runtime.

Many finite element problems require the definition of material parameters.
These can be handled in *SfePy* with material variables which associate the
material parameters with the corresponding region of the mesh.

::

    region_1000 = {
        'name' : 'Omega',
        'select' : 'cells of group 6',
    }

    region_03 = {
        'name' : 'Gamma_Left',
        'select' : 'vertices in (x < 0.00001)',
        'kind' : 'facet',
    }

    region_4 = {
        'name' : 'Gamma_Right',
        'select' : 'vertices in (x > 0.099999)',
        'kind' : 'facet',
    }

Regions assign names to various parts of the finite element mesh. The region
names can later be referred to, for example when specifying portions of the mesh
to apply boundary conditions to. Regions can be specified in a variety of ways,
including by element or by node. Here, Omega is the elemental domain over which
the PDE is solved and Gamma_Left and Gamma_Right define surfaces upon which the
boundary conditions will be applied.

::

    field_1 = {
        'name' : 'temperature',
        'dtype' : 'real',
        'shape' : (1,),
        'region' : 'Omega',
        'approx_order' : 1,
    }

A field is used mainly to define the approximation on a (sub)domain, i.e. to
define the discrete spaces :math:`V_h`, where we seek the solution.

The Poisson equation can be used to compute e.g. a temperature distribution,
so let us call our field ``'temperature'``. On the region ``'Omega'``
it will be approximated using linear finite elements.

A field in a given region defines the finite element approximation.
Several variables can use the same field, see below.

::

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

One field can be used to generate discrete degrees of freedom (DOFs) of
several variables. Here the unknown variable (the temperature) is called
``'t'``, it's associated DOF name is ``'t.0'`` --- this will be referred
to in the Dirichlet boundary section (``ebc``). The corresponding test
variable of the weak formulation is called ``'s'``. Notice that the
``'dual'`` item of a test variable must specify the unknown it
corresponds to.

For each unknown (or state) variable there has to be a test (or virtual)
variable defined, as usual in weak formulation of PDEs.

::

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

Essential (Dirichlet) boundary conditions can be specified as above.

Boundary conditions place restrictions on the finite element formulation and
create a unique solution to the problem. Here, we specify that a temperature of
+2 is applied to the left surface of the mesh and a temperature of -2 is applied
to the right surface.

::

    integral_1 = {
        'name' : 'i',
        'order' : 2,
    }

Integrals specify which numerical scheme to use. Here we are using a 2nd order
quadrature over a 3 dimensional space.

::

    equations = {
        'Temperature' : """dw_laplace.i.Omega( coef.val, s, t ) = 0"""
    }

The equation above directly corresponds to the discrete version of
:eq:`eq_wlaplace_static`, namely:  Find :math:`\bm{t} \in V_h`, such that

.. math::
    \bm{s}^T (\int_{\Omega_h} c\ \bm{G}^T G) \bm{t} = 0, \quad
    \forall \bm{s} \in V_{h0} \;,

where :math:`\nabla u \approx \bm{G} \bm{u}`.

The equations block is the heart of the *SfePy* problem definition file. Here,
we are specifying that the Laplacian of the temperature (in the weak
formulation) is 0, where ``coef.val`` is a material constant. We are using the
``i`` integral defined previously, over the domain specified by the region
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

    solver_0 = {
        'name' : 'ls',
        'kind' : 'ls.scipy_direct',
        'method' : 'auto',
    }

Here, we specify which kind of solver to use for linear equations.

::

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

Even linear problems are solved by a nonlinear solver (KISS rule) - only one
iteration is needed and the final rezidual is obtained for free.

::

    options = {
        'nls' : 'newton',
        'ls' : 'ls',
    }

The solvers to use are specified in the options block. We can define multiple
solvers with different convergence parameters if necessary.

That's it! Now it is possible to proceed as described in
:ref:`invoking_command_line`.

Short syntax of keywords
^^^^^^^^^^^^^^^^^^^^^^^^

The same diffusion equation example as above in **short syntax** reads,
see :download:`examples/diffusion/poisson_short_syntax.py
<../examples/diffusion/poisson_short_syntax.py>`, as follows:

.. literalinclude:: ../examples/diffusion/poisson_short_syntax.py
   :linenos:

.. _sec-interactive-example-linear-elasticity:

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

    \int_{\Omega} D_{ijkl} e_{kl}(\ul{u}) e_{ij}(\ul{v}) + \int_{\Omega}
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
    >>> from sfepy.applications import solve_pde
    >>> from sfepy.postprocess import Viewer

    When in SfePy source directory, try:
    >>> pb, vec, data = solve_pde('examples/diffusion/poisson.py')
    >>> view = Viewer(pb.get_output_name())
    >>> view()

    When in another directory (and SfePy is installed), try:
    >>> from sfepy import data_dir
    >>> pb, vec, data = solve_pde(data_dir + '/examples/diffusion/poisson.py')
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

Define the regions - the whole domain :math:`\Omega`, where the solution
is sought, and :math:`\Gamma_1`, :math:`\Gamma_2`, where the boundary
conditions will be applied. As the domain is rectangular, we first get a
bounding box to get correct bounds for selecting the boundary edges.

.. sourcecode:: ipython

    In [3]: min_x, max_x = domain.get_mesh_bounding_box()[:, 0]
    In [4]: eps = 1e-8 * (max_x - min_x)
    In [5]: omega = domain.create_region('Omega', 'all')
    In [6]: gamma1 = domain.create_region('Gamma1',
       ...:                               'vertices in x < %.10f' % (min_x + eps),
       ...:                               'facet')
    In [7]: gamma2 = domain.create_region('Gamma2',
       ...:                               'vertices in x > %.10f' % (max_x - eps),
       ...:                               'facet')

Next we define the actual finite element approximation using the
:class:`Field` class.

.. sourcecode:: ipython

    In [8]: field = Field.from_args('fu', nm.float64, 'vector', omega,
       ...:                         space='H1', poly_space_base='lagrange',
       ...:                         approx_order=2)

Using the field `fu`, we can define both the unknown variable :math:`\ub` and
the test variable :math:`\vb`.

.. sourcecode:: ipython

    In [9]: u = FieldVariable('u', 'unknown', field)
    In [10]: v = FieldVariable('v', 'test', field, primary_var_name='u')

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

The equations have to be completed by boundary conditions. Let us clamp
the left edge :math:`\Gamma_1`, and shift the right edge
:math:`\Gamma_2` in the :math:`x` direction a bit, depending on the
:math:`y` coordinate.

.. sourcecode:: ipython

   In [20]: from sfepy.fem.conditions import Conditions, EssentialBC
   In [21]: fix_u = EssentialBC('fix_u', gamma1, {'u.all' : 0.0})
   In [22]: def shift_u_fun(ts, coors, bc=None, problem=None, shift=0.0):
      ....:     val = shift * coors[:,1]**2
      ....:     return val
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

And view them.

.. sourcecode:: ipython

    In [32]: view = Viewer('regions.vtk')
    In [33]: view()

You should see this:

.. image:: images/linear_elasticity_regions.png
   :width: 70 %
   :align: center

Finally, we apply the boundary conditions, solve the problem, save and
view the results.

.. sourcecode:: ipython

    In [34]: pb.time_update(ebcs=Conditions([fix_u, shift_u]))
    In [35]: vec = pb.solve()
    In [36]: print nls_status
    In [37]: pb.save_state('linear_elasticity.vtk', vec)
    In [38]: view = Viewer('linear_elasticity.vtk')
    In [39]: view()

This is the resulting image:

.. image:: images/linear_elasticity_solution1.png
   :width: 70 %
   :align: center

The default view is not very fancy. Let us show the displacements by
shifting the mesh. Close the previous window and do:

.. sourcecode:: ipython

    In [56]: view(vector_mode='warp_norm', rel_scaling=2,
       ....:      is_scalar_bar=True, is_wireframe=True)

And the result is:

.. image:: images/linear_elasticity_solution2.png
   :width: 70 %
   :align: center

See the docstring of `view()` and play with its options.

.. _tutorial_interactive_source:

Complete Example as a Script
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The source code: :download:`linear_elasticity.py
<../examples/standalone/interactive/linear_elasticity.py>`. It should be
run from the *SfePy* source directory so that it finds the mesh file.

.. literalinclude:: ../examples/standalone/interactive/linear_elasticity.py
   :linenos:
