User's Guide
============

.. contents:: Table of Contents
   :local:
   :backlinks: top

This manual provides reference documentation to *SfePy* from a user's
perspective.

Running a simulation
--------------------

The following should be run in the top-level directory of the *SfePy* source
tree after compiling the C extension files. See
:ref:`introduction_installation` for full installation instructions info. The
``$`` indicates the command prompt of your terminal.

Basic usage
^^^^^^^^^^^

* ::

    $ ./simple.py examples/diffusion/poisson.py

  * Creates ``cylinder.vtk``

* ::

    $ ./simple.py examples/navier_stokes/stokes.py

  * Creates ``channels_symm944t.vtk``

* ::

    $ ./runTests.py

  * See `Running Tests`_

* ::

    $ ./isfepy

  * Follow the help information printed on startup

Surface extraction
^^^^^^^^^^^^^^^^^^

* ::

    $ ./findSurf.py meshes/quantum/cube.node -

  * Creates ``surf_cube.mesh``

Applications
^^^^^^^^^^^^

* Phononic Materials

  * ::

      $ ./phonon.py -p examples/phononic/band_gaps.py

    * see ``examples/phononic/output/``

* ``schroedinger.py``

  * (order is important below):

    1. ::

        $ ./schroedinger.py --2d --create-mesh

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
                            directory with tests [default: tests]
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

The example problems in the ``examples`` directory can be computed by the script
``simple.py`` which is in the top-level directory of the *SfePy* distribution.
If it is run without arguments, a help message is printed::

    $ ./simple.py
    Usage: simple.py [options] filename_in

    Options:
      --version             show program's version number and exit
      -h, --help            show this help message and exit
      -c "key : value, ...", --conf="key : value, ..."
                            override problem description file items, written as
                            python dictionary without surrouding braces
      -O "key : value, ...", --options="key : value, ..."
                            override options item of problem description, written
                            as python dictionary without surrouding braces
      -o filename           basename of output file(s) [default: <basename of
                            input file>]
      --format=format       output file format, one of: {vtk, h5, mesh} [default:
                            vtk]
      --log=file            log all messages to specified file (existing file will
                            be overwritten!)
      -q, --quiet           do not print any messages to screen
      --save-ebc            save problem state showing EBC (Dirichlet conditions)
      --save-regions        save problem regions as meshes
      --save-regions-as-groups
                            save problem regions in a single mesh but mark them by
                            using different element/node group numbers
      --save-field-meshes   save meshes of problem fields (with extra DOF nodes)
      --solve-not           do not solve (use in connection with --save-*)
      --list=what           list data, what can be one of: {terms}

Additional (stand-alone) examples are in the examples/ directory, e.g.::

    $ python examples/compare_elastic_materials.py

Parametric study example::

    $ ./simple.py examples/diffusion/poisson_parametric_study.py

Common tasks
""""""""""""

* Run a simulation::

    ./simple.py examples/diffusion/poisson.py
    ./simple.py examples/diffusion/poisson.py -o some_results # -> produces some_results.vtk

* Print available terms::

    ./simple.py --list=terms

* Run a simulation and also save Dirichlet boundary conditions::

    ./simple.py --save-ebc examples/diffusion/poisson.py # -> produces an additional .vtk file with BC visualization

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
      --no-offscreen        force no offscreen rendering for --no-show
      --3d                  3d plot mode
      --view=angle,angle[,distance[,focal_point]]
                            camera azimuth, elevation angles, and optionally also
                            distance and focal point coordinates (without []) as
                            in `mlab.view()` [default: if --3d is True: "45,45",
                            else: "0,0"]
      --roll=angle          camera roll angle [default: 0.0]
      --fgcolor=R,G,B       foreground color, that is the color of all text
                            annotation labels (axes, orientation axes, scalar bar
                            labels) [default: 0.0,0.0,0.0]
      --bgcolor=R,G,B       background color [default: 1.0,1.0,1.0]
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
      --wireframe           show wireframe of mesh surface for each data
      --opacity=opacity     global surface and wireframe opacity in [0.0, 1.0]
                            [default: 1.0]
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
      --group-names=name1,...,nameN:...
                            superimpose plots of data in each group
      --subdomains=mat_id_name,threshold_limits,single_color
                            superimpose surfaces of subdomains over each data;
                            example value: mat_id,0,None,True
      --step=step           set the time step [default: 0]
      --anti-aliasing=value
                            value of anti-aliasing [default: mayavi2 default]
      -d 'var_name0,function_name0,par0=val0,par1=val1,...:var_name1,...', --domain-specific='var_name0,function_name0,par0=val0,par1=val1,...:var_name1,...'
                            domain specific drawing functions and configurations

As a simple example, try::

    $ ./simple.py examples/diffusion/poisson.py
    $ ./postproc.py cylinder.vtk

The following window should display:

.. image:: images/postproc_simple.png

The ``-l`` switch lists information contained in a results file, e.g.::

    $ ./postproc.py -l cylinder.vtk
    sfepy: 0: cylinder.vtk
    point scalars
      "node_groups" (354,) range: 0 0 l2_norm_range: 0.0 0.0
        "t" (354,) range: -2.0 2.0 l2_norm_range: 0.0106091 2.0
        cell scalars
          "mat_id" (1348,) range: 6 6 l2_norm_range: 6.0 6.0

.. _sec-problem-description-file:

Problem description file
------------------------

Here we discuss the basic items that users have to specify in their input
files. For complete examples, see the problem description files in the
``examples/`` directory of SfePy.


FE mesh
^^^^^^^

A FE mesh defining a domain geometry can be stored in several formats:

* legacy VTK (``.vtk``)
* custom HDF5 file (``.h5``)
* medit mesh file (``.mesh``)
* tetgen mesh files (``.node``, ``.ele``)
* comsol text mesh file (``.txt``)
* abaqus text mesh file (``.inp``)
* avs-ucd text mesh file (``.inp``)
* hypermesh text mesh file (``.hmascii``)
* hermes3d mesh file (``.mesh3d``)
* nastran text mesh file (``.bdf``)
* gambit neutral text mesh file (``.neu``)
* salome/pythonocc med binary mesh file (``.med``)

Example::

    filename_mesh = 'meshes/3d/cylinder.vtk'

The VTK and HDF5 formats can be used for storing the results. The format
can be selected in options, see :ref:`miscellaneous_options`.

The following geometry elements are supported:

.. image:: images/elements.png

Regions
^^^^^^^

Regions serve to select a certain part of the computational domain (= selection
of nodes and elements of a FE mesh). They are used to define the boundary
conditions, the domains of terms and materials etc.

* Region selection syntax

  * Entity selections

    * ``all``
    * ``nodes of surface``
    * ``nodes of group <integer>``
    * ``nodes of group <str>`` (if mesh format supports reading boundary
      condition nodes)
    * ``nodes in <expr>``
    * ``nodes by <function>``
    * ``node <id>[, <id>, ...]``
    * ``elements of group <integer>``
    * ``elements by <efunction>``
    * ``element <id>[, <id>, ...]`` assumes group 0 (ig = 0)
    * ``element (<ig>, <id>)[, (<ig>, <id>), ...]``
    * ``r.<name of another region>``

  * Notation

    * ``<expr>`` is a logical expression like ``(y <= 0.00001) & (x < 0.11)``
    * ``<function>`` is e.g., ``afunction( x, y, z, otherArgs )``
    * ``<efunction>`` is e.g., ``efunction( domain )``

  * Region operations

    * Node-wise: ``+n``, ``-n``, ``*n`` (union, set difference, intersection)
    * Element-wise: ``+e``, ``-e``, ``*e`` (union, set difference, intersection)

  * Additional specification:

    * 'forbid' : 'group <integer>' - forbid elements of listed groups
    * 'can_cells' : <boolean> - determines whether a region can have cells (volume in 3D)

* Region definition syntax

  * Long syntax: a region is defined by the following Python dictionary
    (denote optional keys)::

        region_<number> = {
            'name' : <name>,
            'select' : <selection>,
            ['forbid'] : group <integer>[, <integer>],
            ['can_cells'] : <boolean>,
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

  * Short syntax::

          regions = {
              <name> : ( <selection>, {[<additional spec.>]} )
          }

    * Example definitions::

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

Fields
^^^^^^

Fields correspond to FE spaces

* Long syntax::

        field_<number> = {
            'name' : <name>,
            'dtype' : <data_type>,
            'shape' : <shape>,
            'region' : <region_name>,
            'approx_order' : <approx_order>
        }

  where
    * <data_type> is a numpy type (float64 or complex128) or
      'real' or 'complex'
    * <shape> is the number of DOFs per node: 1 or (1,) or 'scalar', space
      dimension (2, or (2,) or 3 or (3,)) or 'vector'; it can be other
      positive integer than just 1, 2, or 3
    * <region_name> is the name of region where the field is defined
    * <approx_order> is the FE approximation order, e.g. 0, 1, 2, '1B' (1
      with bubble)

  * Example: scalar P1 elements in 2D on a region Omega::

        field_1 = {
            'name' : 'temperature',
            'dtype' : 'real',
            'shape' : 'scalar',
            'region' : 'Omega',
            'approx_order' : 1
        }

* Short syntax::

          fields = {
              <name> : (<data_type>, <shape>, <region_name>, <approx_order>)
          }

  * Example: scalar P1 elements in 2D on a region Omega::

        fields = {
            'temperature' : ('real', 1, 'Omega', 1),
        }

* The following approximation orders can be used:

  * simplex elements: 1, 2, '1B', '2B'
  * tensor product elements: 0, 1, '1B'

  Optional bubble function enrichment is marked by 'B'.

Variables
^^^^^^^^^

Variables use the FE approximation given by the specified field:

* Long syntax::

        variables_<number> = {
            'name' : <name>,
            'kind' : <kind>,
            'field' : <field_name>,
            ['order' : <order>,]
            ['dual' : <variable_name>,]
            ['history' : <history_size>,]
        }

  where
    * <kind> - 'unknown field', 'test field' or 'parameter field'
    * <order> -  primary variable - order in the global vector of unknowns
    * <history_size> - number of time steps to remember prior to current step

  * Example, long syntax::

        variable_1 = {
            'name' : 't',
            'kind' : 'unknown field',
            'field' : 'temperature',
            'order' : 0, # order in the global vector of unknowns
            'history' : 1,
        }

        variable_2 = {
            'name' : 's',
            'kind' : 'test field',
            'field' : 'temperature',
            'dual' : 't',
        }

* Short syntax::

        variables = {
            <name> : (<kind>, <field_name>, <spec.>, [<history>])
        }

  where

  * <spec> - in case of: primary variable - order in the global vector of unknowns, dual variable - name of primary variable


  * Example, short syntax::

        variables = {
            't' : ('unknown field', 'temperature', 0, 1),
            's' : ('test field', 'temperature', 't'),
        }

.. _ug_integrals:

Integrals
^^^^^^^^^

Define the integral type and quadrature rule. This keyword is optional.

* Long syntax::

        integral_<number> = {
            'name' : <name>,
            'kind' : <kind>,
            'quadrature' : <rule>
        }

  where

    * <name> - the integral name - it has to begin with 'i'!
    * <kind> - volume 'v' or surface 's' integral
    * <rule> - <family>_o<order>_d<dimension>, available quadratures are in sfe/fem/quadratures.py - it is still preliminary and incomplete

  * Example, long syntax::

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

* Short syntax::

        integrals = {
            <name> : (<kind>, <rule>)
        }

  * Example, short syntax::

        import numpy as nm
        N = 2
        integrals = {
            'i1' : ('v', 'gauss_o2_d3'),
            'i2' : ('v', 'custom', zip(nm.linspace( 1e-10, 0.5, N ),
                                       nm.linspace( 1e-10, 0.5, N )),
                    [1./N] * N),
        }

Boundary conditions
^^^^^^^^^^^^^^^^^^^

The boundary conditions apply in a given region given by its name, and,
optionally, in selected times. The times can be given either using a
list of tuples `(t0, t1)` making the condition active for `t0 <= t <
t1`, or by a name of a function taking the time argument and returning
True or False depending on whether the condition is active at the given
time or not.

* Dirichlet (essential) boundary conditions, long syntax::

    ebc_<number> = {
        'name' : <name>,
        'region' : <region_name>,
        ['times' : <times_specification>,]
        'dofs' : {<dof_specification> : <value>[,
                  <dof_specification> : <value>, ...]}
    }

  * Example::

        ebc_1 = {
            'name' : 'ZeroSurface',
            'region' : 'Surface',
            'times' : [(0.5, 1.0), (2.3, 5)],
            'dofs' : {'u.all' : 0.0, 'phi.all' : 0.0},
        }

* Dirichlet (essential) boundary conditions, short syntax::

    ebcs = {
        <name> : (<region_name>, [<times_specification>,]
                  {<dof_specification> : <value>[,
                   <dof_specification> : <value>, ...]},...)
    }

  * Example::

        ebcs = {
            'u1' : ('Left', {'u.all' : 0.0}),
            'u2' : ('Right', [(0.0, 1.0)], {'u.0' : 0.1}),
            'phi' : ('Surface', {'phi.all' : 0.0}),
        }

Initial conditions
^^^^^^^^^^^^^^^^^^

Initial conditions are applied prior to the boundary conditions - no special
care must be used for the boundary dofs.

* Long syntax::

    ic_<number> = {
        'name' : <name>,
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
        <name> : (<region_name>, {<dof_specification> : <value>[,
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

* Example, long syntax::

    material_10 = {
        'name' : 'm',
        'values' : {
            # This gets tiled to all physical QPs (constant function)
            'val' : [0.0, -1.0, 0.0],
            # This does not - '.' denotes a special value, e.g. a flag.
            '.val0' : [0.0, 0.1, 0.0],
        },
    }

    material_3 = {
      'name' : 'm2',
      'function' : 'some_function',
    }

    def some_function(ts, coor, region, ig, mode=None):
        out = {}
        if mode == 'qp':
            # <array of shape (coor.shape[0], n_row, n_col)>
            out['val'] = nm.ones((coor.shape[0], 1, 1), dtype=nm.float64)
        else: # special mode
            out['val0'] = True

* Example, short syntax::

    material = {
        'm' : ({'val' : [0.0, -1.0, 0.0]},),
        'm2' : 'some_function',
        'm3' : (None, 'some_function'), # Same as the above line.
    }


Equations and Terms
^^^^^^^^^^^^^^^^^^^

Equations can be built by combining terms listed in :ref:`term_table`.

Examples
""""""""

* Laplace equation, named integral::

    equations = {
        'Temperature' : """dw_laplace.i1.Omega( coef.val, s, t ) = 0"""
    }

* Laplace equation, simplified integral given by order::

    equations = {
        'Temperature' : """dw_laplace.2.Omega( coef.val, s, t ) = 0"""
    }

* Laplace equation, automatic integration order (not implemented yet!)::

    equations = {
        'Temperature' : """dw_laplace.a.Omega( coef.val, s, t ) = 0"""
    }

* Navier-Stokes equations::

    equations = {
        'balance' :
        """+ dw_div_grad.i2.Omega( fluid.viscosity, v, u )
           + dw_convect.i2.Omega( v, u )
           - dw_stokes.i1.Omega( v, p ) = 0""",
        'incompressibility' :
        """dw_stokes.i1.Omega( u, q ) = 0""",
    }

Configuring Solvers
^^^^^^^^^^^^^^^^^^^

In SfePy, a non-linear solver has to be specified even when solving a linear
problem. The linear problem is/should be then solved in one iteration of the
nonlinear solver.

* Linear solver, long syntax::

    solver_0 = {
        'name' : 'ls',
        'kind' : 'ls.umfpack',
    }

* Nonlinear solver, long syntax::

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


* Solvers, short syntax::

    solvers = {
        'ls' : ('ls.scipy_direct', {}),
        'newton' : ('nls.newton',
                    {'i_max'   : 1,
                     'problem' : 'nonlinear'}),
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
function keyword. See the examples below.

Defining material parameters
""""""""""""""""""""""""""""

The functions for defining material parameters can work in two modes,
distinguished by the `mode` argument. The two modes are 'qp' and 'special'. The
first mode is used for usual functions that define parameters in quadrature
points (hence 'qp'), while the second one can be used for special values like
various flags.

The shape and type of data returned in the 'special' mode can be arbitrary
(depending on the term used). On the other hand, in the 'qp' mode all the data
have to be numpy float64 arrays with shape `(n_coor, n_row, n_col)`, where
`n_coor` is the number of quadrature points given by the `coors` argument,
`n_coor = coors.shape[0]`, and `(n_row, n_col)` is the shape of a material
parameter in each quadrature point. For example, for scalar parameters, the
shape is `(n_coor, 1, 1)`.

Examples
""""""""

See ``examples/diffusion/poisson_functions.py`` for a complete problem
description file demonstrating how to use different kinds of functions.

- functions for defining regions::

    def get_circle(coors, domain=None):
        r = nm.sqrt(coors[:,0]**2.0 + coors[:,1]**2.0)
        return nm.where(r < 0.2)[0]

    functions = {
        'get_circle' : (get_circle,),
    }

- functions for defining boundary conditions::

    def get_p_edge(ts, coors, bc=None, problem=None):
        if bc.name == 'p_left':
            return nm.sin(nm.pi * coors[:,1])
        else:
            return nm.cos(nm.pi * coors[:,1])

    functions = {
        'get_p_edge' : (get_p_edge,),
    }

    ebcs = {
        'p' : ('Gamma', {'p.0' : 'get_p_edge'}),
    }

  The values can be given by a function of time, coordinates and
  possibly other data, for example::

    ebcs = {
        'f1' : ('Gamma1', {'u.0' : 'get_ebc_x'}),
        'f2' : ('Gamma2', {'u.all' : 'get_ebc_all'}),
    }

    def get_ebc_x(coors, amplitude):
        z = coors[:, 2]
        val = amplitude * nm.sin(z * 2.0 * nm.pi)
        return val

    def get_ebc_all(ts, coors):
        x, y, z = coors[:, 0], coors[:, 1], coors[:, 2]
        val = ts.step * nm.r_[x, y, z]
        return val

    functions = {
        'get_ebc_x' : (lambda ts, coors, bc, problem, **kwargs:
                       get_ebc_x(coors, 5.0),),
        'get_ebc_all' : (lambda ts, coors, bc, problem, **kwargs:
                         get_ebc_all(ts, coors),),
    }

  Note that when setting more than one component as in `get_ebc_all()`
  above, the function should return a single one-dimensional vector with
  all values of the first component, then of the second one
  etc. concatenated together.

- function for defining usual material parameters::

    def get_pars(ts, coors, mode=None, region=None, ig=None):
        if mode == 'qp':
            val = coors[:,0]
            val.shape = (coors.shape[0], 1, 1)

            return {'x_coor' : val}

    functions = {
        'get_pars' : (get_pars,),
    }

- function for defining special material parameters, with an extra argument::

    def get_pars_special(ts, coors, mode=None, region=None, ig=None,
                         extra_arg=None):
        if mode == 'special':
            if extra_arg == 'hello!':
                ic = 0
            else:
                ic = 1
            return {('x_%s' % ic) : coors[:,ic]}

    functions = {
        'get_pars1' : (lambda ts, coors, mode=None, region=None, ig=None:
                       get_pars_special(ts, coors, mode, region, ig,
                                        extra_arg='hello!'),),
    }

    # Just another way of adding a function, besides 'functions' keyword.
    function_1 = {
        'name' : 'get_pars2',
        'function' : lambda ts, coors,mode=None,  region=None, ig=None:
            get_pars_special(ts, coors, mode, region, ig, extra_arg='hi!'),
    }

- function combining both kinds of material parameters::

    def get_pars_both(ts, coors, mode=None, region=None, ig=None):
        out = {}

        if mode == 'special':

            out['flag'] = coors.max() > 1.0

        elif mode == 'qp':

            val = coors[:,1]
            val.shape = (coors.shape[0], 1, 1)

            out['y_coor'] = val

        return out

    functions = {
        'get_pars_both' : (get_pars_both,),
    }

- function for setting values of a parameter variable::

    variable_1 = {
        'name' : 'p',
        'kind' : 'parameter field',
        'field' : 'temperature',
        'like' : None,
        'special' : {'setter' : 'get_load_variable'},
    }

    def get_load_variable(ts, coors, region=None):
        y = coors[:,1]
        val = 5e5 * y
        return val

    functions = {
        'get_load_variable' : (get_load_variable,)
    }

.. _miscellaneous_options:

Miscellaneous
^^^^^^^^^^^^^
The options can be used to select solvers, output file format, output
directory, to register functions to be called at various phases of the
solution (the `hooks`), and for other settings.

* Additional options (including solver selection)::

    options = {
        # string, output directory
        'output_dir'        : 'output/<output_dir>',

        # 'vtk' or 'h5', output file (results) format
        'output_format'     : 'h5',

        # string, nonlinear solver name
        'nls' : 'newton',

        # string, linear solver name
        'ls' : 'ls',

        # string, time stepping solver name
        'ts' : 'ts',

        # int, number of time steps when results should be saved (spaced
        # regularly from 0 to n_step), or -1 for all time steps
        'save_steps' : -1,

        # string, a function to be called after each time step
        'step_hook'  : '<step_hook_function>',

        # string, a function to be called after each time step, used to
        # update the results to be saved
        'post_process_hook' : '<post_process_hook_function>',

        # string, as above, at the end of simulation
        'post_process_hook_final' : '<post_process_hook_final_function>',

        # string, a function to generate probe instances
        'gen_probes'        : '<gen_probes_function>',

        # string, a function to probe data
        'probe_hook'        : '<probe_hook_function>',

        # string, a function to modify problem definition parameters
        'parametric_hook' : '<parametric_hook_function>',
    }

  * ``post_process_hook`` enables computing derived quantities, like
    stress or strain, from the primary unknown variables. See the
    examples in ``examples/large_deformation/`` directory.
  * ``parametric_hook`` makes it possible to run parametric studies by
    modifying the problem description programmatically. See
    ``examples/diffusion/poisson_parametric_study.py`` for an example.
  * ``output_dir`` redirects output files to specified directory

Building Equations in SfePy
---------------------------

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

Available Solvers
-----------------

This Section describes solvers available in SfePy from user's
perspective. There internal/external solvers include linear, nonlinear,
eigenvalue, optimization and time stepping solvers.

Nonlinear Solvers
^^^^^^^^^^^^^^^^^

Almost every problem, even linear, is solved in SfePy using a nonlinear
solver that calls a linear solver in each iteration. This approach
unifies treatment of linear and non-linear problems, and simplifies
application of Dirichlet (essential) boundary conditions, as the linear
system computes not a solution, but a solution increment, i.e., it
always has zero boundary conditions.

The following solvers are available:

- 'nls.newton': Newton solver with backtracking line-search - this is
  the default solver, that is used for almost all examples.
- 'nls.oseen': Oseen problem solver tailored for stabilized
  Navier-Stokes equations (see :ref:`navier_stokes-stabilized_navier_stokes`).
- 'nls.scipy_broyden_like': interface to Broyden and Anderson solvers
  from scipy.optimize.
- 'nls.semismooth_newton': Semismooth Newton method for contact/friction
  problems.

Linear solvers
^^^^^^^^^^^^^^

A good linear solver is key to solving efficiently stationary as well as
transient PDEs with implicit time-stepping. The following solvers are
available:

- 'ls.scipy_direct': direct solver from SciPy - this is the default
  solver for all examples. It is strongly recommended to install umfpack
  and its SciPy wrappers to get good performance.
- 'ls.umfpack': alias to 'ls.scipy_direct'.
- 'ls.scipy_iterative': Interface to SciPy iterative solvers.
- 'ls.pyamg': Interface to PyAMG solvers.
- 'ls.petsc': Interface to Krylov subspace solvers of PETSc.
- 'ls.petsc_parallel': Interface to Krylov subspace solvers of PETSc
  able to run in parallel by storing the system to disk and running a
  separate script via `mpiexec`.
- 'ls.schur_complement': Schur complement problem solver.
