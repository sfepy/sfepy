"""
Common code for basic electronic structure examples.

Notes
-----

The same code should work also with a 3D (box) mesh, but a very fine mesh would
be required. Also in the 2D case, finer mesh and/or higher approximation order
means higher accuracy.

Try changing C, F and L parameters in square.geo and regenerate the mesh using
gmsh::

  $ gmsh -2 -format mesh meshes/quantum/square.geo -o meshes/quantum/square.mesh
  $ ./script/convert_mesh.py meshes/quantum/square.mesh meshes/quantum/aux.vtk
  $ ./script/convert_mesh.py meshes/quantum/aux.vtk meshes/quantum/square.mesh

The ``script/convert_mesh.py`` calls make the mesh 2D, as gmsh does not save
planar medit meshes.

Also try changing approximation order ('approx_order') of the field below, as
well as the integral order (should be two times the approximation order).
"""
from __future__ import absolute_import
from sfepy import data_dir

def common(fun_v, n_eigs=5, tau=0.0):
    filename_mesh = data_dir + '/meshes/quantum/square.mesh'

    options = {
        'save_eig_vectors' : None,
        'n_eigs' : n_eigs,
        'eigen_solver' : 'eigen1',
    }

    region_1000 = {
        'name' : 'Omega',
        'select' : 'all',
    }

    region_2 = {
        'name' : 'Surface',
        'select' : 'vertices of surface',
        'kind' : 'facet',
    }

    functions = {
        'fun_v' : (fun_v,),
    }

    material_1 = {
        'name' : 'm',

        'values' : {
            'val' : 0.5,
        },
    }

    material_2 = {
        'name' : 'mat_v',

        'function' : 'fun_v',
    }

    field_0 = {
        'name' : 'field_Psi',
        'dtype' : 'real',
        'shape' : 'scalar',
        'region' : 'Omega',
        'approx_order' : 2,
    }

    integral_1 = {
        'name' : 'i',
        'order' : 4,
    }

    variable_1 = {
        'name' : 'Psi',
        'kind' : 'unknown field',
        'field' : 'field_Psi',
        'order' : 0,
    }
    variable_2 = {
        'name' : 'v',
        'kind' : 'test field',
        'field' : 'field_Psi',
        'dual' : 'Psi',
    }
    variable_3 = {
        'name' : 'V',
        'kind' : 'parameter field',
        'field' : 'field_Psi',
        'like' : 'Psi',
    }

    ebc_1 = {
        'name' : 'ZeroSurface',
        'region' : 'Surface',
        'dofs' : {'Psi.0' : 0.0},
    }

    equations = {
        'lhs' : """  dw_laplace.i.Omega( m.val, v, Psi )
                   + dw_volume_dot.i.Omega( mat_v.V, v, Psi )""",
        'rhs' : """dw_volume_dot.i.Omega( v, Psi )""",
    }

    solver_2 = {
        'name' : 'eigen1',
        'kind' : 'eig.pysparse',

        'tau' : tau,
        'eps_a' : 1e-10,
        'i_max' : 150,
        'method' : 'qmrs',
        'verbosity' : 0,
        'strategy' : 1,
    }

    return locals()
