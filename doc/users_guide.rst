User's Guide
============

This manual provides reference documentation to *SfePy* from a user's perspective.

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
  
      $ ./eigen.py -p examples/phononic/band_gaps.py
      
    * see ``examples/phononic/output/``

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
* medit mesh file (``.mesh``)
* tetgen mesh files (``.node``, ``.ele``)
* comsol text mesh file (``.txt``)
* custom HDF5 file (``.h5``)

Example::

    filename_mesh = 'meshes/3d/cylinder.vtk'

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
    * 'canCells' : <boolean> - determines whether a region can have cells (volume in 3D) 

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
            'dim' : (<dofs>,1),
            'domain' : <region_name>,
            'bases' : {<subregion_name> : <bases>}
	    ['dtype' : <dtype>]
    	}


  where
    * <dofs> - number of DOFs per node
    * <dtype> - 'real' or 'complex' values
    * <bases> - approximation on subdomains, e.g.  {'Omega_1' : 3_4_P1, 'Omega_2' : '3_8_Q1'} 
	
  * Example: P1 elements in 2D on a whole domain Omega::

        field_1 = {
            'name' : 'temperature',
            'dim' : (1,1),
            'domain' : 'Omega',
            'bases' : {'Omega' : '2_3_P1'}
        }
	
* Short syntax::

     	  fields = {
	      <name> : ((<dofs>,1), <dtype>, <region_name>, {<subregion_name> : <bases>})
	  }

* The following elements/approximations can be used:

  * 2D: 2_3_P1, 2_3_P1B, 2_3_P2, 2_3_P2B, 2_4_Q0, 2_4_Q1
  * 3D: 3_4_P0, 3_4_P1, 3_4_P1B, 3_4_P2, 3_4_P2B, 3_8_Q0, 3_8_Q1

The letter 'P' indicates a polynomial space on the simplex geometry, while 'Q'
on the tensor product geometry. The numbers behind the letters indicate
approximation order. Optional bubble function enrichment is marked by 'B'.

.. image:: images/elements.png

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
    	}

  where
    * <kind> - 'unknown field', 'test field' or 'parameter field'
    * <order> -  primary variable - order in the global vector of unknowns

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

* Short syntax::

    	variables = {
            <name> : (<kind>, <field_name>, <spec.>)
        }

  where

  * <spec> - in case of: primary variable - order in the global vector of unknowns, dual variable - name of primary variable


  * Example, short syntax::

        variables = {
            't' : ('unknown field', 'temperature', 0),
            's' : ('test field', 'temperature', 't'),
        }

Integrals
^^^^^^^^^

Define the integral type and quadrature rule. 

* Long syntax::
    
        integral_<number> = {
            'name' : <name>,
	    'kind' : <kind>,
            'quadrature' : <rule>
    	}

  where

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

* Dirichlet (essential) boundary conditions, long syntax::

    ebc_<number> = {
        'name' : <name>,
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
        <name> : (<region_name>, {<dof_specification> : <value>[,
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
            # <array of shape (coor.shape[0], n_row, n_col)>
            out['val'] = nm.ones((coor.shape[0], 1, 1), dtype=nm.float64)
        else: # special mode
            out['val0'] = True

* Example, short syntax::

    material = {
        'm' : ('SomeRegion', {'val' : [0.0, -1.0, 0.0]}),
    }


Equations and Terms
^^^^^^^^^^^^^^^^^^^

Equations can be built by combining terms listed in :ref:`term_table`.

Examples
""""""""

* Laplace equation::

    equations = {
        'Temperature' : """dw_laplace.i1.Omega( coef.val, s, t ) = 0"""
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

- functions for defining boundary conditions (`get_p_edge()`) and regions
  (`get_circle()`)::

    def get_p_edge(ts, coors, bc=None):
        if bc.name == 'p_left':
            return nm.sin(nm.pi * coors[:,1])
        else:
            return nm.cos(nm.pi * coors[:,1])

    def get_circle(coors, domain=None):
        r = nm.sqrt(coors[:,0]**2.0 + coors[:,1]**2.0)
        return nm.where(r < 0.2)[0]

    functions = {
        'get_p_edge' : (get_p_edge,),
        'get_circle' : (get_circle,),
    }

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

