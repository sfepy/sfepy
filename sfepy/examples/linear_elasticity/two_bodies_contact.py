r"""
Contact of two elastic bodies with a penalty function for enforcing the contact
constraints.

Find :math:`\ul{u}` such that:

.. math::
    \int_{\Omega} D_{ijkl}\ e_{ij}(\ul{v}) e_{kl}(\ul{u})
    + \int_{\Gamma_{c}} \varepsilon_N \langle g_N(\ul{u}) \rangle \ul{n} \ul{v}
    = 0
    \;, \quad \forall \ul{v} \;,

where :math:`\varepsilon_N \langle g_N(\ul{u}) \rangle` is the penalty
function, :math:`\varepsilon_N` is the normal penalty parameter, :math:`\langle
g_N(\ul{u}) \rangle` are the Macaulay's brackets of the gap function
:math:`g_N(\ul{u})` and

.. math::
    D_{ijkl} = \mu (\delta_{ik} \delta_{jl}+\delta_{il} \delta_{jk}) +
    \lambda \ \delta_{ij} \delta_{kl}
    \;.

Usage examples::

  sfepy-run sfepy/examples/linear_elasticity/two_bodies_contact.py --save-regions-as-groups --save-ebc-nodes

  sfepy-view two_bodies.mesh.vtk -f u:wu:f2:p0 1:vw:p0 gap:p1 -2

  python3 sfepy/scripts/plot_logs.py log.txt

  sfepy-view two_bodies.mesh_ebc_nodes.vtk -2
  sfepy-view two_bodies.mesh_regions.vtk -2
"""
from sfepy.mechanics.matcoefs import stiffness_from_youngpoisson
from sfepy.discrete.fem.meshio import UserMeshIO

import numpy as nm

dim = 2

if dim == 2:
    dims0 = [1.0, 0.5]
    shape0 = [4, 4]
    centre0 = [0, -0.25]

    dims1 = [1.0, 0.5]
    shape1 = [3, 3]
    centre1 = [0, 0.25]

    shift1 = [0.0, -0.1]

else:
    dims0 = [1.0, 1.0, 0.5]
    shape0 = [2, 2, 2]
    centre0 = [0, 0, -0.25]

    dims1 = [1.0, 1.0, 0.5]
    shape1 = [2, 2, 2]
    centre1 = [0, 0, 0.25]

    shift1 = [0.0, 0.0, -0.1]

def get_bbox(dims, centre, eps=0.0):
    dims = nm.asarray(dims)
    centre = nm.asarray(centre)

    bbox = nm.r_[[centre - (0.5 - eps) * dims], [centre + (0.5 - eps) * dims]]
    return bbox

def gen_two_bodies(dims0, shape0, centre0, dims1, shape1, centre1, shift1):
    from sfepy.discrete.fem import Mesh
    from sfepy.mesh.mesh_generators import gen_block_mesh

    m0 = gen_block_mesh(dims0, shape0, centre0)
    m1 = gen_block_mesh(dims1, shape1, centre1)

    coors = nm.concatenate((m0.coors, m1.coors + shift1), axis=0)

    desc = m0.descs[0]
    c0 = m0.get_conn(desc)
    c1 = m1.get_conn(desc)
    conn = nm.concatenate((c0, c1 + m0.n_nod), axis=0)

    ngroups = nm.zeros(coors.shape[0], dtype=nm.int32)
    ngroups[m0.n_nod:] = 1

    mat_id = nm.zeros(conn.shape[0], dtype=nm.int32)
    mat_id[m0.n_el:] = 1

    name = 'two_bodies.mesh'

    mesh = Mesh.from_data(name, coors, ngroups, [conn], [mat_id], m0.descs)

    return mesh

def mesh_hook(mesh, mode):
    if mode == 'read':
        return gen_two_bodies(dims0, shape0, centre0,
                              dims1, shape1, centre1, shift1)

    elif mode == 'write':
        pass

def post_process(out, pb, state, extend=False):
    from sfepy.base.base import Struct
    from sfepy.discrete.fem import extend_cell_data

    ev = pb.evaluate
    gap = ev('dw_contact.i.Contact(contact.epss, v, u)',
             mode='el_avg', term_mode='gap')
    gap = extend_cell_data(gap, pb.domain, 'Contact', val=0.0, is_surface=True)
    out['gap'] = Struct(name='output_data',
                        mode='cell', data=gap, dofs=None)

    return out

filename_mesh = UserMeshIO(mesh_hook)

options = {
    'nls' : 'newton',
    'ls' : 'ls',
    'post_process_hook' : 'post_process',
}

fields = {
    'displacement': ('real', dim, 'Omega', 1),
}

materials = {
    'solid' : ({'D': stiffness_from_youngpoisson(dim,
                                                 young=1.0, poisson=0.3)},),
    'contact' : ({'.epss' : 1e1},),
}

variables = {
    'u' : ('unknown field', 'displacement', 0),
    'v' : ('test field', 'displacement', 'u'),
}

bbox0 = get_bbox(dims0, centre0, eps=1e-5)
bbox1 = get_bbox(dims1, nm.asarray(centre1) + nm.asarray(shift1), eps=1e-5)

if dim == 2:
    regions = {
        'Omega' : 'all',
        'Omega0' : 'cells of group 0',
        'Omega1' : 'cells of group 1',
        'Bottom' : ('vertices in (y < %f)' % bbox0[0, 1], 'facet'),
        'Top' : ('vertices in (y > %f)' % bbox1[1, 1], 'facet'),
        'Contact0' : ('(vertices in (y > %f) *v r.Omega0)' % bbox0[1, 1],
                      'facet'),
        'Contact1' : ('(vertices in (y < %f) *v r.Omega1)' % bbox1[0, 1],
                      'facet'),
        'Contact' : ('r.Contact0 +s r.Contact1', 'facet')
    }

else:
    regions = {
        'Omega' : 'all',
        'Omega0' : 'cells of group 0',
        'Omega1' : 'cells of group 1',
        'Bottom' : ('vertices in (z < %f)' % bbox0[0, 2], 'facet'),
        'Top' : ('vertices in (z > %f)' % bbox1[1, 2], 'facet'),
        'Contact0' : ('(vertices in (z > %f) *v r.Omega0)' % bbox0[1, 2],
                      'facet'),
        'Contact1' : ('(vertices in (z < %f) *v r.Omega1)' % bbox1[0, 2],
                      'facet'),
        'Contact' : ('r.Contact0 +s r.Contact1', 'facet')
    }

ebcs = {
    'fixb' : ('Bottom', {'u.all' : 0.0}),
    'fixt' : ('Top', {'u.all' : 0.0}),
}

integrals = {
    'i' : 10,
}

equations = {
    'elasticity' :
    """dw_lin_elastic.2.Omega(solid.D, v, u)
     + dw_contact.i.Contact(contact.epss, v, u)
     = 0""",
}

solvers = {
    'ls' : ('ls.scipy_direct', {}),
    'newton' : ('nls.newton', {
            'i_max' : 5,
            'eps_a' : 1e-6,
            'eps_r' : 1.0,
            'macheps' : 1e-16,
            # Linear system error < (eps_a * lin_red).
            'lin_red' : 1e-2,
            'ls_red' : 0.1,
            'ls_red_warp' : 0.001,
            'ls_on' : 100.1,
            'ls_min' : 1e-5,
            'check' : 0,
            'delta' : 1e-8,
            # 'log' : {'text' : 'log.txt', 'plot' : None},
    })
}
