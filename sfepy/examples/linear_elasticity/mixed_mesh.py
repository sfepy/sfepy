r"""
A linear elastic beam loaded with a continuous force. The FE mesh consists
of hexehedral and tetrahedral elements.

The displacement at the beam end is compared to the reference
solution calculated on the homogeneous hexahedral mesh.

Running the simulation::

    sfepy-run sfepy/examples/linear_elasticity/mixed_mesh.py

Viewing the results::

    sfepy-view beam_h* -f u:s0:wu:e:p0 u:s1:wu:e:p0 --camera-position="1.2,-0.6,0.1,0.4,0.1,-0.1,-0.2,0.1,1"
"""
from __future__ import absolute_import
from sfepy import data_dir
from sfepy.mechanics.matcoefs import stiffness_from_youngpoisson
from sfepy.base.base import Struct
import numpy as nm


def reference_solution(pb):
    from sfepy.discrete import Problem

    conf = pb.conf.copy()
    conf.filename_mesh = data_dir + '/meshes/3d/beam_h7.mesh'
    for k in list(conf.regions.keys()):
        if '_t' in k:
            del conf.regions[k]

    for k in list(conf.fields.keys()):
        if '_t' in k:
            del conf.fields[k]

    for k in list(conf.variables.keys()):
        if '_t' in k:
            del conf.variables[k]

    conf.equations = {
        'balance_of_forces':
            """dw_lin_elastic.i.Omega(solid.D, v, u)
             = dw_surface_ltr.i.Top(force.val, v)"""
    }

    rpb = Problem.from_conf(conf)
    rpb.set_output_dir(pb.output_dir)

    return rpb.solve().get_state_parts()['u']


def get_force(ts, coors, mode=None, **kwargs):
    if mode == 'qp':
        F = 1e3

        val = nm.zeros_like(coors)[..., None]
        val[:, 2, 0] = -coors[:, 0] / 0.7 * F

        return {'val': val}


def post_proces(out, pb, state, extend=False):
    c1 = pb.domain.regions['Omega'].cells
    c2 = pb.domain.regions['Omega_t'].cells
    ed = pb.domain.regions['Edge'].vertices

    S = nm.zeros((pb.domain.cmesh.n_el, 1, 6, 1), dtype=nm.float64)

    S[c1] = pb.evaluate('ev_cauchy_stress.i.Omega(solid.D, u)',
                        mode='el_avg')
    S[c2] = pb.evaluate('ev_cauchy_stress.i.Omega_t(solid.D, u_t)',
                        mode='el_avg')

    out['stress'] = Struct(name='out', data=S, mode='cell')

    u_nd = out['u'].data[ed[0], :]
    ref_u_nd = reference_solution(pb).reshape(-1, 3)[ed[0], :]

    du = nm.linalg.norm((ref_u_nd - u_nd) / nm.linalg.norm(ref_u_nd))
    print(f'Relative difference with respect to the reference solution: {du}')

    return out


filename_mesh = data_dir + '/meshes/3d/beam_h5t12.mesh'

options = {
    'post_process_hook': 'post_proces',
}

regions = {
    'Omega_': 'all',
    'Omega': ('cells of group 1', 'cell', None, {'vertices_from': 'Omega_'}),
    'Omega_t': ('cells of group 2', 'cell', None, {'vertices_from': 'Omega_'}),
    'Left': ('vertices in (x < 0.01)', 'facet'),
    'Right': ('vertices in (x > 0.69)', 'facet'),
    'Top_': ('vertices in (z > 0.09)', 'facet'),
    'Bottom': ('vertices in (z < 0.01)', 'facet'),
    'Edge': ('r.Bottom *v r.Right', 'vertex'),
    'Top': ('r.Top_ *f r.Omega', 'facet'),
    'Top_t': ('r.Top_ *f r.Omega_t', 'facet'),
}

functions = {
    'get_force' : (get_force,),
}

materials = {
    'solid': ({'D': stiffness_from_youngpoisson(dim=3, young=1e6, poisson=0.3)},),
    'force': 'get_force',
    'vforce': ({'.val' : [0.0, 0.0, -7]},),
}

fields = {
    'displacement': ('real', 'vector', 'Omega', 1),
    'displacement_t': ('real', 'vector', 'Omega_t', 1),
}

integrals = {
    'i': 2,
}

variables = {
    'u': ('unknown field', 'displacement', 0),
    'v': ('test field', 'displacement', 'u'),
    'u_t': ('unknown field', 'displacement_t', 0),
    'v_t': ('test field', 'displacement_t', 'u_t'),
}

ebcs = {
    'Fixed': ('Left', {'u.all' : 0.0}),
}

equations = {
    'balance_of_forces':
    """dw_lin_elastic.i.Omega(solid.D, v, u)
     + dw_lin_elastic.i.Omega_t(solid.D, v_t, u_t)
     = dw_surface_ltr.i.Top(force.val, v)
     + dw_surface_ltr.i.Top_t(force.val, v_t)"""
}

solvers = {
    'ls': ('ls.auto_direct', {}),
    'newton': ('nls.newton', {
        'i_max': 1,
        'eps_a': 1e-6,
    }),
}
