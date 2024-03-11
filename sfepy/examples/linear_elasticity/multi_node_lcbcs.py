# -*- coding: utf-8 -*-
r"""
Linear elasticity with multi node linear combination constraints.

Find :math:`\ul{u}` such that:

.. math::
    \int_{\Omega} D_{ijkl}\ e_{ij}(\ul{v}) e_{kl}(\ul{u})
    = 0
    \;, \quad \forall \ul{v} \;,

where

.. math::
    D_{ijkl} = \mu (\delta_{ik} \delta_{jl}+\delta_{il} \delta_{jk}) +
    \lambda \ \delta_{ij} \delta_{kl}
    \;.

The hanging nodes in ``'Dependent_nodes'`` region are bounded
to the nodes in ``'Independent_nodes'`` region using the ``lcbcs`
(``multi_node_combination``) conditions

View the results using::

  sfepy-view hanging_nodes.vtk -f u:wu:e -2
"""
import numpy as nm
from sfepy import data_dir
from sfepy.base.base import output
from sfepy.mechanics.matcoefs import stiffness_from_lame


def node_map_fun(coors0, coors1):
    """
    Map hanging nodes (coors0) and coincident edges (edge nodes in coors1)
    on a line.
    """

    ldir = coors1[-1] - coors1[0]
    ldir = ldir / nm.linalg.norm(ldir)
    if ldir[0] > 1e-9:
        par0 = coors0[:, 0] / ldir[0]
        par1 = coors1[:, 0] / ldir[0]
    else:
        par0 = coors0[:, 1] / ldir[1]
        par1 = coors1[:, 1] / ldir[1]

    idxs1 = nm.argsort(par1)
    nmap = []
    for par in par0:
        edge = nm.where(nm.logical_and(par1[:-1] < par, par < par1[1:]))[0][0]
        de = par1[edge + 1] - par1[edge]
        nmap.append((idxs1[edge], idxs1[edge + 1],
                     1 - (par - par1[edge]) / de,
                     1 - (par1[edge + 1] - par) / de))

    return (nm.arange(par0.shape[0]),
            nm.array([k[:2] for k in nmap]),
            nm.array([k[2:] for k in nmap]))


def check_lcs(out, pb, state, extend=None):
    displ = state.get_state_parts()['u'].reshape((-1, 2))
    coors = pb.domain.mesh.coors

    hnodes_map = {12: [2, 6], 14: [2, 6], 16: [6, 10], 18: [6, 10]}

    ok = True
    for nd, ends in hnodes_map.items():
        dy = coors[ends[1], 1] - coors[ends[0], 1]
        t = (coors[nd, 1] - coors[ends[0], 1]) / dy
        u_exp = displ[ends[0]] * (1 - t) + displ[ends[1]] * t
        ok = ok and (nm.linalg.norm(displ[nd] - u_exp) < 1e-9)

    result = {True: 'passed', False: 'failed'}[ok]
    output(f'multi node linear condition test {result}!')

    # tests_lcbcs hook:
    nls_status = pb.get_solver().status.nls_status
    if hasattr(nls_status, 'conditions'):
        nls_status.conditions.append(int(not ok))

    return out


filename_mesh = data_dir + '/meshes/2d/special/hanging_nodes.mesh'

options = {
    'post_process_hook': check_lcs,
}

regions = {
    'Omega': 'all',
    'Omega1': 'cells of group 1',
    'Omega2': 'cells of group 2',
    'Inside_nodes': ('vertices in (x > 0.199) & (x < 0.201)', 'vertex'),
    'Independent_nodes': ('r.Inside_nodes *v r.Omega1', 'vertex'),
    'Dependent_nodes': ('r.Inside_nodes -v r.Independent_nodes', 'vertex'),  
    'Left': ('vertices in (x < 0.001)', 'facet'),
    'Right': ('vertices in (x > 0.299)', 'facet'),
}

materials = {
    'solid': ({'D': stiffness_from_lame(dim=2, lam=1e1, mu=1e0)},),
}

fields = {
    'displacement': ('real', 'vector', 'Omega', 1),
}

integrals = {
    'i': 1,
}

variables = {
    'u': ('unknown field', 'displacement'),
    'v': ('test field', 'displacement', 'u'),
}

ebcs = {
    'Fixed': ('Left', {'u.all': 0.0}),
    'Displaced': ('Right', {'u.0': 0.05, 'u.1': -0.05}),
}

functions = {
    'node_map_fun': (node_map_fun,),
}

lcbcs = {
    'hanging': (['Dependent_nodes', 'Independent_nodes'], {'u.all': 'u.all'},
                'node_map_fun', 'multi_node_combination', None),
}

equations = {
    'balance_of_forces' :
    """dw_lin_elastic.i.Omega(solid.D, v, u) = 0""",
}

solvers = {
    'ls': ('ls.auto_direct', {}),
    'newton': ('nls.newton', {
        'i_max': 1,
        'eps_a': 1e-10,
    }),
}
