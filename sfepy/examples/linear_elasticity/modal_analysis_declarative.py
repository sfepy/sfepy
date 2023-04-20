"""
Modal analysis of a 3D cylinder.

The first six modes are the rigid body modes because no boundary
conditions are applied.

Running the simulation::

  sfepy-run sfepy/examples/linear_elasticity/modal_analysis_declarative.py

The eigenvalues are saved to cylinder_eigs.txt and the eigenvectros to
cylinder.vtk. View the results using::

  sfepy-view cylinder.vtk -f u008:wu008:f0.005
"""
import numpy as nm
from sfepy.base.base import output
from sfepy.mechanics.matcoefs import stiffness_from_youngpoisson
from sfepy import data_dir


def report_eigs(pb, evp):
    eigs = evp.eigs
    n_rbm = evp.pb.conf.n_rbm

    output('%d eigenvalues converged (%d ignored as rigid body modes)' %
           (len(eigs), n_rbm))

    eigs = eigs[n_rbm:]

    omegas = nm.sqrt(eigs)
    freqs = omegas / (2 * nm.pi)

    output('number |         eigenvalue |  angular frequency '
           '|          frequency')
    for ii, eig in enumerate(eigs):
        output('%6d | %17.12e | %17.12e | %17.12e'
               % (ii + 1, eig, omegas[ii], freqs[ii]))


def define(n_eigs=5, approx_order=1):
    filename_mesh, dim = data_dir + '/meshes/3d/cylinder.mesh', 3

    n_rbm = dim * (dim + 1) // 2  # number of rigid body modes

    options = {
        'n_eigs': n_eigs + n_rbm,
        'eigs_only': False,
        'post_process_hook_final': 'report_eigs',
        'evps': 'eig',
    }

    regions = {
        'Omega': 'all',
    }

    materials = {
        'm': ({
            'D': stiffness_from_youngpoisson(dim, 210e+9, 0.3),
            'rho': 7850.0,
        },),
    }

    fields = {
        'displacement': ('real', 'vector', 'Omega', approx_order),
    }

    variables = {
        'u': ('unknown field', 'displacement'),
        'v': ('test field', 'displacement', 'u'),
    }

    integrals = {
        'i': 2 * approx_order,
    }

    equations = {
        'lhs': 'dw_lin_elastic.i.Omega(m.D, v, u)',
        'rhs': 'dw_dot.i.Omega(m.rho, v, u)',
    }

    ebcs = {}

    solvers = {
        'eig': ('eig.scipy', {
            'method': 'eigsh',
            'tol': 1e-3,
            'maxiter': 1000,
        }),
    }

    return locals()
