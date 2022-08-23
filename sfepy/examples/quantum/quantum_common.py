"""
Common code for basic electronic structure examples.

It covers only simple single electron problems, e.g. well, oscillator, hydrogen
atom and boron atom with 1 electron - see the corresponding files in this
directory, where potentials (:func:`fun_v()`) as well as exact solutions
(:func:`get_exact()`) for those problems are defined.

Notes
-----

The same code should work also with a 3D (box) mesh, but a very fine mesh would
be required. Also in the 2D case, finer mesh and/or higher approximation order
means higher accuracy.

Try changing C, F and L parameters in ``meshes/quantum/square.geo`` and
regenerate the mesh using gmsh::

  gmsh -2 -format mesh meshes/quantum/square.geo -o meshes/quantum/square.mesh
  ./script/convert_mesh.py -2 meshes/quantum/square.mesh meshes/quantum/square.mesh

The ``script/convert_mesh.py`` call makes the mesh planar, as gmsh saves 2D
medit meshes including the zero z coordinates.

Also try changing approximation order ('approx_order') of the field below.

Usage Examples
--------------

The following examples are available and can be run using::

  sfepy-run sfepy/examples/quantum/boron.py
  sfepy-run sfepy/examples/quantum/hydrogen.py
  sfepy-run sfepy/examples/quantum/oscillator.py
  sfepy-run sfepy/examples/quantum/well.py
"""
from __future__ import absolute_import
from sfepy.base.base import output
from sfepy import data_dir

def common(fun_v, get_exact=None, n_eigs=5, tau=0.0):

    def report_eigs(pb, evp):
        from numpy import NaN

        bounding_box = pb.domain.mesh.get_bounding_box()
        box_size = bounding_box[1][0] - bounding_box[0][0]
        output('box_size: %f' % box_size)
        output('eigenvalues:')

        if get_exact is not None:
            eeigs = get_exact(n_eigs, box_size, pb.domain.shape.dim)

            output('n      exact         FEM      error')
            for ie, eig in enumerate(evp.eigs):
                if ie < len(eeigs):
                    exact = eeigs[ie]
                    err = 100*abs((exact - eig)/exact)
                else:
                    exact = NaN
                    err = NaN
                output('%d:  %.8f   %.8f  %7.4f%%' % (ie, exact, eig, err))

        else:
            output('n       FEM')
            for ie, eig in enumerate(evp.eigs):
                output('%d:  %.8f' % (ie, eig))

    filename_mesh = data_dir + '/meshes/quantum/square.mesh'

    options = {
        'n_eigs' : n_eigs,
        'eigs_only' : False,
        'post_process_hook_final' : 'report_eigs',

        'evps' : 'eig',
    }

    regions = {
        'Omega' : 'all',
        'Surface' : ('vertices of surface', 'facet'),
    }

    materials = {
        'm' : ({'val' : 0.5},),
        'mat_v' : 'fun_v',
    }

    functions = {
        'fun_v' : (fun_v,),
    }

    approx_order = 2
    fields = {
        'field_Psi' : ('real', 'scalar', 'Omega', approx_order),
    }

    variables = {
        'Psi' : ('unknown field', 'field_Psi', 0),
        'v' : ('test field', 'field_Psi', 'Psi'),
    }

    ebcs = {
        'ZeroSurface' : ('Surface', {'Psi.0' : 0.0}),
    }

    integrals = {
        'i' : 2 * approx_order,
    }

    equations = {
        'lhs' : """dw_laplace.i.Omega(m.val, v, Psi)
                 + dw_dot.i.Omega(mat_v.V, v, Psi)""",
        'rhs' : """dw_dot.i.Omega(v, Psi)""",
    }

    solvers = {
        'eig' : ('eig.scipy', {
            'method' : 'eigsh',
            'tol' : 1e-10,
            'maxiter' : 150,

            # Compute the eigenvalues near tau using the shift-invert mode.
            'which' : 'LM',
            'sigma' : tau,
        }),
    }

    return locals()
