Tutorial
========

*SfePy* can be used in two basic ways:
  #. a black-box partial differential equation (PDE) solver,
  #. a Python package to build custom applications involving solving PDEs by the
     finite element (FE) method.

This tutorial focuses on the first way, and introduces the basic concepts and
nomenclature used in the following parts of the documentation.

Notes on solving PDEs by the Finite element method
--------------------------------------------------

Equations in *SfePy* are built using terms, which correspond directly to the
integral forms of weak formulation of a problem to be solved. As an example, let
us consider the Laplace equation in time interval :math:`t \in [0, t_{\rm
final}]`: 

.. math::
   :label: eq_laplace

    \pdiff{T}{t} + c \Delta T = 0 \mbox{ in }\Omega,\quad T(t) = \bar{T}(t)
    \mbox{ on } \Gamma \;.

The weak formulation of :eq:`eq_laplace` is: Find :math:`T \in V`, such that 

.. math:: 
   :label: eq_wlaplace

    \int_{\Omega} s \pdiff{T}{t} + \int_{\Omega} c\ \nabla T : \nabla s = 0,
    \quad \forall s \in V_0 \;,

where we assume no fluxes over :math:`\partial \Omega \setminus \Gamma`. In the
syntax used in *SfePy* input files, this can be written as::

    dw_mass_scalar.i1.Omega( s, dT/dt ) + dw_laplace.i1.Omega( coef, s, T) = 0

which directly corresponds to the discrete version of :eq:`eq_wlaplace`: Find
:math:`\bm{T} \in V_h`, such that 

.. math::

    \bm{s}^T (\int_{\Omega_h} \bm{\phi}^T \bm{\phi}) \pdiff{\bm{T}}{t} +
    \bm{s}^T (\int_{\Omega_h} c\ \bm{G}^T \bm{G}) \bm{T} = 0, \quad \forall
    \bm{s} \in V_{h0} \;,

where :math:`u \approx \bm{\phi} \bm{u}`, :math:`\nabla u \approx \bm{G}
\bm{u}` for :math:`u \in \{s, T\}`. The integrals over the discrete domain
:math:`\Omega_h` are approximated by a numerical quadrature, that is named
:math:`\verb|i1|` in our case.

Term call syntax
^^^^^^^^^^^^^^^^

In general, the syntax of a term call in *SfePy* is::

    <term_name>.<i>.<r>( <arg1>, <arg2>, ... )

where ``<i>`` denotes an integral name (i.e. a name of numerical quadrature to
use) and ``<r>`` marks a region (domain of the integral). In the following,
``<virtual>`` corresponds to a test function, ``<state>`` to a unknown function
and ``<parameter>`` to a known function arguments.

Running a simulation
--------------------

The following should be run in the top-level directory of the *SfePy* source
tree after compiling the C extension files. See
:ref:`introduction_installation` for full installation instructions info. The
``$`` indicates the command prompt of your terminal.

Basic usage
^^^^^^^^^^^

* ::

    $ ./simple.py input/poisson.py

  * Creates ``simple.vtk``

* ::

    $ ./simple.py input/pfdpm_permeability.py

  * Creates ``perf_symm944t.vtk``

* ::

    $ ./runTests.py

  * See `Running Tests`_

* ::

    $ ./isfepy

  * Follow the help information printed on startup 

Surface extraction
^^^^^^^^^^^^^^^^^^

* ::

    $ ./findSurf.py database/t.1.node

  * Creates ``surf_t.1.mesh``

Applications
^^^^^^^^^^^^

* Phononic Materials

  * ::
  
      $ ./eigen.py -p input/phono/band_gaps.py
      
    * see ``input/phono/output/``

* ``schroedinger.py`` 

  * (order is important below):

    1. :: 

        $ ./schroedinger.py --2d --mesh

    2. :: 

        $ ./schroedinger.py --2d --hydrogen

    3. ::
    
        $ ./postproc.py mesh.vtk

Stand-Alone Examples
^^^^^^^^^^^^^^^^^^^^

* :: 

    $ python examples/rs_correctors.py

* :: 

    $ python examples/compare_elastic_materials.py

* ::

    $ python examples/live_plot.py

Running Tests
^^^^^^^^^^^^^

The tests are run by the ``runTests.py`` script::

    $ ./runTests.py -h
    Usage: runTests.py [options] [test_filename[ test_filename ...]]

    Options:
      --version             show program's version number and exit
      -h, --help            show this help message and exit
      --print-doc           print the docstring of this file (howto write new
                            tests)
      -d directory, --dir=directory
                            directory with tests [default: <top_dir>/tests]
      -o directory, --output=directory
                            directory for storing test results and temporary files
                            [default: output-tests]
      --debug               raise silenced exceptions to see what was wrong
      --filter-none         do not filter any messages
      --filter-less         filter output (suppress all except test messages)
      --filter-more         filter output (suppress all except test result
                            messages)

Common tasks
""""""""""""

* Run all tests, filter output; result files related to the tests can be found
  in output-tests directory::

    ./runTests.py
    ./runTests.py --filter-more
    ./runTests.py --filter-less

* Run a particular test file, filter output::

    ./runTests.py tests/test_input_le.py # Test if linear elasticity input file works.

* Debug a failing test::

    ./runTests.py tests/test_input_le.py --debug

Computations and examples
^^^^^^^^^^^^^^^^^^^^^^^^^

The example problems in the ``input`` directory can be computed by the script
``simple.py`` which is in the top-level directory of the *SfePy* distribution.
If it is run without arguments, a help message is printed::

    $ ./simple.py
    Usage: simple.py [options] filename_in

    Options:
      --version             show program's version number and exit
      -h, --help            show this help message and exit
      -o filename           basename of output file(s) [default: <basename of
                            input file>]
      --format=format       output file format, one of: {vtk, h5, mesh} [default:
                            vtk]
      --save-ebc            save problem state showing EBC (Dirichlet conditions)
      --save-regions        save problem regions as meshes
      --save-field-meshes   save meshes of problem fields (with extra DOF nodes)
      --save-region-field-meshes
                            save meshes of regions of problem fields (with extra
                            DOF nodes)
      --solve-not           do not solve (use in connection with --save-*)
      --list=what           list data, what can be one of: {terms}

Additional (stand-alone) examples are in the examples/ directory, e.g.::

    $ python examples/compare_elastic_materials.py

Parametric study example::

    $ ./simple.py input/poisson_parametric_study.py

Common tasks
""""""""""""

* Run a simulation::

    ./simple.py input/poisson.py
    ./simple.py input/poisson.py -o some_results # -> produces some_results.vtk

* Print available terms::

    ./simple.py --list=terms

* Run a simulation and also save Dirichlet boundary conditions::

    ./simple.py --save-ebc input/poisson.py # -> produces an additional .vtk file with BC visualization

Visualization of results
------------------------

The ``postproc.py`` script can be used for quick postprocessing and
visualization of the *SfePy* results. It requires mayavi2 installed on your
system. Running ``postproc.py`` without arguments produces::

    $ ./postproc.py                        
    Usage: postproc.py [options] filename                                        

    This is a script for quick Mayavi-based visualizations of finite element
    computations results.                                                   

    Examples
    --------
      The examples assume that runTests.py has been run successfully and the
      resulting data files are present.                                     

      - view data in output-tests/test_navier_stokes.vtk

        $ python postproc.py output-tests/test_navier_stokes.vtk
        $ python postproc.py output-tests/test_navier_stokes.vtk --3d

      - create animation (forces offscreen rendering) from
        output-tests/test_time_poisson.*.vtk              
                                                          
        $ python postproc.py output-tests/test_time_poisson.*.vtk -a mov

      - create animation (forces offscreen rendering) from
        output-tests/test_hyperelastic.*.vtk              

        The range specification for the displacements 'u' is required, as
        output-tests/test_hyperelastic.00.vtk contains only zero         
        displacements which leads to invisible glyph size.               
                                                                         
        $ python postproc.py output-tests/test_hyperelastic.*.vtk                          --ranges=u,0,0.02 -a mov 

      - same as above, but slower frame rate

        $ python postproc.py output-tests/test_hyperelastic.*.vtk                          --ranges=u,0,0.02 -a mov --ffmpeg-options="-r 2 -sameq"



    Options:
      --version             show program's version number and exit
      -h, --help            show this help message and exit       
      -l, --list-ranges     do not plot, only list names and ranges of all data
      -n, --no-show         do not call mlab.show()                            
      --3d                  3d plot mode                                       
      --view=angle,angle    camera view angles [default: if --3d is True: "45,45",
                            else: "0,0"]                                          
      --roll=angle          camera roll angle [default: 0.0]                      
      --layout=layout       layout for multi-field plots, one of: rowcol, colrow, 
                            row, col [default: rowcol]                            
      --scalar-mode=mode    mode for plotting scalars with --3d, one of:          
                            cut_plane, iso_surface, both [default: iso_surface]   
      --vector-mode=mode    mode for plotting vectors, one of: arrows, norm,      
                            arrows_norm, warp_norm [default: arrows_norm]         
      -s scale, --scale-glyphs=scale                                              
                            relative scaling of glyphs (vector field              
                            visualization) [default: 0.05]                        
      --clamping            glyph clamping mode                                   
      --ranges=name1,min1,max1:name2,min2,max2:...                                
                            force data ranges [default: automatic from data]      
      -b, --scalar-bar      show scalar bar for each data                         
      --rel-text-width=width                                                      
                            relative text annotation width [default: 0.02]        
      -w, --watch           watch the results file for changes (single file mode  
                            only)
      -o filename, --output=filename
                            view image file name [default: 'view.png']
      --output-dir=directory
                            output directory for saving view images; ignored when
                            -o option is given, as the directory part of the
                            filename is taken instead [default: '.']
      -a <ffmpeg-supported format>, --animation=<ffmpeg-supported format>
                            if set to a ffmpeg-supported format (e.g. mov, avi,
                            mpg), ffmpeg is installed and results of multiple time
                            steps are given, an animation is created in the same
                            directory as the view images
      --ffmpeg-options="<ffmpeg options>"
                            ffmpeg animation encoding options (enclose in "")
                            [default: -r 10 -sameq]
      -r resolution, --resolution=resolution
                            image resolution in NxN format [default: shorter axis:
                            600; depends on layout: for rowcol it is 800x600]
      --all                 draw all data (normally, node_groups and mat_id are
                            omitted)
      --only-names=list of names
                            draw only named data
      --step=step           set the time step [default: 0]
      --anti-aliasing=value
                            value of anti-aliasing [default: mayavi2 default]

As a simple example, try::

    $ ./simple.py input/poisson.py
    $ ./postproc.py simple.vtk

The following window should display:

.. image:: images/postproc_simple.png

The ``-l`` switch lists information contained in a results file, e.g.::

    $ ./postproc.py -l simple.vtk
    sfepy: 0: simple.vtk
    point scalars
      "node_groups" (354,) range: 0 0 l2_norm_range: 0.0 0.0
        "t" (354,) range: -2.0 2.0 l2_norm_range: 0.0106091 2.0
        cell scalars
          "mat_id" (1348,) range: 6 6 l2_norm_range: 6.0 6.0

Problem description file
------------------------

Here we discuss the basic items that users have to specify in their input
files. For complete examples, see the problem description files in the
``input/`` directory of SfePy.


FE mesh
^^^^^^^

A FE mesh defining a domain geometry can be stored in several formats:

* legacy VTK (``.vtk``)
* medit mesh file (``.mesh``)
* tetgen mesh files (``.node``, ``.ele``)
* comsol text mesh file (``.txt``)
* custom HDF5 file (``.h5``)

Example::

    filename_mesh = 'database/a_mesh.vtk'

Regions
^^^^^^^

Regions serve to select a certain part of the computational domain (= selection
of nodes and elements of a FE mesh). They are used to define the boundary
conditions, the domains of terms and materials etc.

* Region selection syntax

  * Entity selections

    * ``all``
    * ``nodes of surface``
    * ``nodes in <expr>``
    * ``nodes by <function>``
    * ``node <id>``
    * ``elements of group <integer>``
    * ``elements by <efunction>``
    * ``r.<name of another region>``

  * Notation

    * ``<expr>`` is a logical expression like ``(y <= 0.00001) & (x < 0.11)``
    * ``<function>`` is e.g., ``afunction( x, y, z, otherArgs )``
    * ``<efunction>`` is e.g., ``efunction( domain )``

  * Region operations

    * Node-wise: ``+n``, ``-n``, ``*n`` (union, set difference, intersection)
    * Element-wise: ``+e``, ``-e``, ``*e`` (union, set difference, intersection)

* Region definition syntax

  * Long syntax: a region is defined by the following Python dictionary
    (denote optional keys)::

        region_<number> = {
            'name' : <string>,
            'select' : <selection string>,
            ['forbid'] : group <integer>[, <integer>]* # forbid elements of listed groups
            ['can_cells'] : <boolean> # determines whether a region can have cells (volume in 3D)
        }

    * Example definitions::

            region_20 = {
                'name' : 'Left',
                'select' : 'nodes in (x < -0.499)'
            }
            region_21 = {
                'name' : 'Right',
                'select' : 'nodes in (x > 0.499)'
            }
            region_31 = {
                'name' : 'Gamma1',
                'select' : """(elements of group 1 *n elements of group 4)
                              +n
                              (elements of group 2 *n elements of group 4)
                              +n
                              ((r.Left +n r.Right) *n elements of group 4)
                           """,
                'forbid' : 'group 1 2'
            }

  * Example definitions, short syntax::

        regions = {
            'Left' : ('nodes in (x < -0.499)', {}),
            'Right' : ('nodes in (x > 0.499)', {}),
            'Gamma1' : ("""(elements of group 1 *n elements of group 4)
                           +n
                           (elements of group 2 *n elements of group 4)
                           +n
                           ((r.Left +n r.Right) *n elements of group 4)""",
                         {'forbid' : 'group 1 2'}),
        }

Fields, Variables and Integrals
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Fields correspond to FE spaces

  * Example: P1 elements in 2D on a whole domain Omega::

        field_1 = {
            'name' : 'temperature',
            'dim' : (1,1),
            'domain' : 'Omega',
            'bases' : {'Omega' : '2_3_P1'}
        }

* Variables use the FE approximation given by the specified field:
          
  * Example, long syntax::

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

  * Example, short syntax::

        variables = {
            't' : ('unknown field', 'temperature', 0),
            's' : ('test field', 'temperature', 't'),
        }

* Integrals (quadrature rules)::

    integral_1 = {
        'name' : 'i1',
        'kind' : 'v',
        'quadrature' : 'gauss_o2_d2', # <quadrature name>
    }

    import numpy as nm
    N = 2
    integral_2 = {
        'name' : 'i2',
        'kind' : 'v',
        'quadrature' : 'custom', # <quadrature name>
        'vals'    : zip(nm.linspace( 1e-10, 0.5, N ),
                        nm.linspace( 1e-10, 0.5, N )),
        'weights' : [1./N] * N,
    }

  * Available quadratures are in sfe/fem/quadratures.py - it is still
    preliminary and incomplete
  * Naming convention: ``<family>_o<order>_d<dimension>``

Boundary conditions
^^^^^^^^^^^^^^^^^^^

* Dirichlet (essential) boundary conditions, long syntax::

    ebc_<number> = {
        'name' : <string>,
        'region' : <region_name>,
        'dofs' : {<dof_specification> : <value>[,
                  <dof_specification> : <value>, ...]}
    }

  * Example::

        ebc_1 = {
            'name' : 'ZeroSurface',
            'region' : 'Surface',
            'dofs' : {'u.all' : 0.0, 'phi.all' : 0.0},
        }

* Dirichlet (essential) boundary conditions, short syntax::

    ebcs = {
        'name' : (<region_name>, {<dof_specification> : <value>[,
                                  <dof_specification> : <value>, ...]},...)
    }

  * Example::

        ebcs = {
            'u1' : ('Left', {'u.all' : 0.0}),
            'u2' : ('Right', {'u.0' : 0.1}),
            'phi' : ('Surface', {'phi.all' : 0.0}),
        }

Initial conditions
^^^^^^^^^^^^^^^^^^

Initial conditions are applied prior to the boundary conditions - no special
care must be used for the boundary dofs.

* Long syntax::

    ic_<number> = {
        'name' : <string>,
        'region' : <region_name>,
        'dofs' : {<dof_specification> : <value>[,
                  <dof_specification> : <value>, ...]}
    }

  * Example::

        ic_1 = {
            'name' : 'ic',
            'region' : 'Omega',
            'dofs' : {'T.0' : 5.0},
        }

* Short syntax::

    ics = {
        'name' : (<region_name>, {<dof_specification> : <value>[,
                                  <dof_specification> : <value>, ...]},...)
    }

  * Example::

        ics = {
            'ic' : ('Omega', {'T.0' : 5.0}),
        }

Materials
^^^^^^^^^

Materials are used to define constitutive parameters (e.g. stiffness,
permeability, or viscosity), and other non-field arguments of terms (e.g. known
traction or volume forces). Depending on a particular term, the parameters can
be constants, functions defined over FE mesh nodes, functions defined in the
elements, etc.

* Example::

    material_10 = {
        'name' : 'm',
        'region' : 'SomeRegion',
        'values' : {
            # This gets tiled to all physical QPs (constant function)
            'val' : [0.0, -1.0, 0.0],
            # This does not - '.' denotes a special value, e.g. a flag.
            '.val0' : [0.0, 0.1, 0.0],
        },
    }

    material_3 = {
      'name' : 'm',
      'region' : 'SomeRegion',
      'function' : 'some_function',
    }

    def some_function(ts, coor, region, ig, mode=None):
        out = {}
        if mode == 'qp':
            out['val'] = <array of shape (coor.shape[0], n_row, n_col)>
        else: # special mode
            out['val0'] = True

Equations and Terms
^^^^^^^^^^^^^^^^^^^

Equations can be built by combining terms listed in ``sfepy_manual.pdf``.

Examples
""""""""

* Laplace equation::

    equations = {
        'Temperature' : """dw_laplace.i1.Omega( coef.val, s, t ) = 0"""
    }

* Navier-Stokes equations::

    equations = {
        'balance' :
        """+ dw_div_grad.i2.Omega( fluid.viscosity, v, w )
           + dw_convect.i2.Omega( v, w )
           - dw_grad.i1.Omega( v, r ) = 0""",
        'incompressibility' :
        """dw_div.i1.Omega( q, w ) = 0""",
    }

Configuring Solvers
^^^^^^^^^^^^^^^^^^^

In SfePy, a non-linear solver has to be specified even when solving a linear
problem. The linear problem is/should be then solved in one iteration of the
nonlinear solver.

* Linear solver::

    solver_0 = {
        'name' : 'ls',
        'kind' : 'ls.umfpack',
    }

* Nonlinear solver::

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

* Solver selection::

    options = {
        'nls' : 'newton',
        'ls' : 'ls',
    }

Functions
^^^^^^^^^

Functions are a way of customizing *SfePy* behavior. They make it possible to
define material properties, boundary conditions, parametric sweeps, and other
items in an arbitrary manner. Functions are normal Python functions declared in
the Problem Definition file, so they can invoke the full power of Python. In
order for *SfePy* to make use of the functions, they must be declared using the
function keyword. See below for examples::

    def get_pars(ts, coors, mode=None, region=None, ig=None, extra_arg=None):
        if mode == 'special':
            if extra_arg == 'hello!':
                ic = 0
            else:
                ic = 1
            return {('x_%s' % ic) : coors[:,ic]}

    def get_p_edge(ts, coors, bc=None):
        if bc.name == 'p_left':
            return nm.sin(nm.pi * coors[:,1])
        else:
            return nm.cos(nm.pi * coors[:,1])

    def get_circle(coors, domain=None):
        r = nm.sqrt(coors[:,0]**2.0 + coors[:,1]**2.0)
        return nm.where(r < 0.2)[0]

    functions = {
        'get_pars1' : (lambda ts, coors, mode=None, region=None, ig=None:
                       get_pars(ts, coors, mode, region, ig, extra_arg='hello!'),),
        'get_p_edge' : (get_p_edge,),
        'get_circle' : (get_circle,),
    }

    # Just another way of adding a function, besides 'functions' keyword.
    function_1 = {
        'name' : 'get_pars2',
        'function' : lambda ts, coors,mode=None,  region=None, ig=None:
            get_pars(ts, coors, mode, region, ig, extra_arg='hi!'),
    }

Miscellaneous
^^^^^^^^^^^^^

* Number of elements assembled in one term function call::

    fe = {
        'chunk_size' : 1000
    }

* Additional options (including solver selection)::

    options = {
        'parametric_hook' : '<a_function>',
        'output_dir' : 'output/<output_dir>',
        'nls' : 'newton',
        'ls' : 'ls',
    }

  * ``parametric_hook`` makes it possible to run parametric studies by
    modifying the problem description programmatically. See
    ``input/poisson_parametric_study.py`` for an example.
  * ``output_dir`` redirects output files to specified directory

