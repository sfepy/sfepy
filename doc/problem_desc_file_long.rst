
Problem Description File - Long Syntax
--------------------------------------

Historically, the keywords exist in two flavors:

- long syntax is the original one - it is longer to type, but the
  individual fields are named

- short syntax was added later to offer brevity.

Region Definition Syntax
^^^^^^^^^^^^^^^^^^^^^^^^

Region, long syntax::

      region_<number> = {
          'name' : <name>,
          'select' : <selection>,
          ['kind'] : <region kind>,
          ['parent'] : <parent region>,
      }

**Example definitions**::

          region_0 = {
              'name' : 'Omega',
              'select' : 'all',
          }
          region_21 = {
              'name' : 'Right',
              'select' : 'vertices in (x > 0.99)',
              'kind' : 'facet',
          }
          region_31 = {
              'name' : 'Gamma1',
              'select' : """(cells of group 1 *v cells of group 2)
                            +v r.Right""",
              'kind' : 'facet',
              'parent' : 'Omega',
          }


Fields
^^^^^^

Fields, long syntax::

       field_<number> = {
            'name' : <name>,
            'dtype' : <data_type>,
            'shape' : <shape>,
            'region' : <region_name>,
            'approx_order' : <approx_order>,
            ['space' : <space>,]
            ['poly_space_basis' : <poly_space_basis>,]
        }

see :ref:`User's Guide-Fields` for meaning of <data_type>, <shape>,
<region_name>, <approx_order>, <space> and <>poly_space_basis>.

**Example**: scalar P1 elements in 2D on a region Omega::

        field_1 = {
            'name' : 'temperature',
            'dtype' : 'real',
            'shape' : 'scalar',
            'region' : 'Omega',
            'approx_order' : 1
        }


Variables
^^^^^^^^^

Variables, long syntax::

        variables_<number> = {
            'name' : <name>,
            'kind' : <kind>,
            'field' : <field_name>,
            ['order' : <order>,]
            ['dual' : <variable_name>,]
            ['history' : <history>,]
        }

where
   * <kind> - 'unknown field', 'test field' or 'parameter field'
   * <order> -  primary variable - order in the global vector of unknowns
   * <history> - number of time steps to remember prior to current step

**Example**::

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


Integrals
^^^^^^^^^

Integrals, long syntax::

        integral_<number> = {
            'name' : <name>,
            'order' : <order>,
        }

where
    * <name> - the integral name - it has to begin with 'i'!
    * <order> - the order of polynomials to integrate, or 'custom' for
      integrals with explicitly given values and weights

**Example**::

        integral_1 = {
            'name' : 'i1',
            'order' : 2,
        }

        import numpy as nm
        N = 2
        integral_2 = {
            'name' : 'i2',
            'order' : 'custom',
            'vals'    : zip(nm.linspace( 1e-10, 0.5, N ),
                            nm.linspace( 1e-10, 0.5, N )),
            'weights' : [1./N] * N,
        }

Essential Boundary Conditions and Constraints
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

See :ref:`User's Guide-EssentialBC` for details.

Dirichlet Boundary Conditions
"""""""""""""""""""""""""""""

Dirichlet (essential) boundary conditions, long syntax::

    ebc_<number> = {
        'name' : <name>,
        'region' : <region_name>,
        ['times' : <times_specification>,]
        'dofs' : {<dof_specification> : <value>[,
                  <dof_specification> : <value>, ...]}
    }

**Example**::

        ebc_1 = {
            'name' : 'ZeroSurface',
            'region' : 'Surface',
            'times' : [(0.5, 1.0), (2.3, 5)],
            'dofs' : {'u.all' : 0.0, 'phi.all' : 0.0},
        }

Periodic Boundary Conditions
""""""""""""""""""""""""""""

Periodic boundary conditions, long syntax::

    epbc_<number> = {
        'name' : <name>,
        'region' : (<region1_name>, <region2_name>),
        ['times' : <times_specification>,]
        'dofs' : {<dof_specification> : <dof_specification>[,
                  <dof_specification> : <dof_specification>, ...]},
        'match' : <match_function_name>,
    }

**Example**::

        epbc_1 = {
            'name' : 'up1',
            'region' : ('Left', 'Right'),
            'dofs' : {'u.all' : 'u.all', 'p.0' : 'p.0'},
            'match' : 'match_y_line',
        }

Linear Combination Boundary Conditions
""""""""""""""""""""""""""""""""""""""

Linear combination boundary conditions, long syntax::

    lcbc_<number> = {
        'name' : <name>,
        'region' : (<region1_name>, <region2_name>) | <region1_name>,
        ['times' : <times_specification>,]
        'dofs' : {<dof_specification> : <dof_specification> | None[, ...]},
        ['dof_map_fun' : <dof_map_function_name> | None,]
        'kind' : <lcbc_kind>,
        [<kind_specific_options>]
    }

**Example**::

        lcbc_1 = {
            'name' : 'rigid',
            'region' : 'Y2',
            'dofs' : {'u.all' : None},
            'kind' : 'rigid',
        }

Initial Conditions
^^^^^^^^^^^^^^^^^^

Initial conditions, long syntax::

    ic_<number> = {
        'name' : <name>,
        'region' : <region_name>,
        'dofs' : {<dof_specification> : <value>[,
                  <dof_specification> : <value>, ...]}
    }

**Example**::

        ic_1 = {
            'name' : 'ic',
            'region' : 'Omega',
            'dofs' : {'T.0' : 5.0},
        }

Materials
^^^^^^^^^

**Example**::

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
      'function' : 'get_pars',
    }

    def get_pars(ts, coors, mode=None, **kwargs):
        out = {}
        if mode == 'qp':
            # <array of shape (coors.shape[0], n_row, n_col)>
            out['val'] = nm.ones((coors.shape[0], 1, 1), dtype=nm.float64)
        else: # special mode
            out['val0'] = True

        return out

Configuring Solvers
^^^^^^^^^^^^^^^^^^^

Linear solver::

    solver_0 = {
        'name' : 'ls',
        'kind' : 'ls.scipy_direct',
    }

Nonlinear solver::

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
        'is_linear' : False,
    }
