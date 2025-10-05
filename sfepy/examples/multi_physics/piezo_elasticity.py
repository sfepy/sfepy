r"""
Piezo-elasticity problem - linear elastic material with piezoelectric
effects.

Find :math:`\ul{u}`, :math:`\phi` such that:

.. math::
    - \omega^2 \int_{Y} \rho\ \ul{v} \cdot \ul{u}
    + \int_{Y} D_{ijkl}\ e_{ij}(\ul{v}) e_{kl}(\ul{u})
    - \int_{Y_2} g_{kij}\ e_{ij}(\ul{v}) \nabla_k \phi
    = 0
    \;, \quad \forall \ul{v} \;,

    \int_{Y_2} g_{kij}\ e_{ij}(\ul{u}) \nabla_k \psi
    + \int_{Y} K_{ij} \nabla_i \psi \nabla_j \phi
    = 0
    \;, \quad \forall \psi \;,

where

.. math::
    D_{ijkl} = \mu (\delta_{ik} \delta_{jl}+\delta_{il} \delta_{jk}) +
    \lambda \ \delta_{ij} \delta_{kl}
    \;.
"""
import os
import numpy as nm

from sfepy import data_dir
from sfepy.discrete.fem import MeshIO
from sfepy.mechanics.matcoefs import stiffness_from_lame

def post_process(out, pb, state, extend=False):
    """
    Calculate and output the strain and stresses for the given state.
    """
    from sfepy.base.base import Struct
    from sfepy.discrete.fem import extend_cell_data

    ev = pb.evaluate
    strain = ev('ev_cauchy_strain.i.Y(u)', mode='el_avg')
    stress = ev('ev_cauchy_stress.i.Y(inclusion.D, u)', mode='el_avg')

    piezo = -ev('ev_piezo_stress.i.Y2(inclusion.coupling, phi)',
                mode='el_avg')
    piezo = extend_cell_data(piezo, pb.domain, 'Y2', val=0.0)

    piezo_strain = ev('ev_piezo_strain.i.Y(inclusion.coupling, u)',
                      mode='el_avg')

    out['cauchy_strain'] = Struct(name='output_data', mode='cell',
                                  data=strain, dofs=None)
    out['elastic_stress'] = Struct(name='output_data', mode='cell',
                                   data=stress, dofs=None)
    out['piezo_stress'] = Struct(name='output_data', mode='cell',
                                 data=piezo, dofs=None)
    out['piezo_strain'] = Struct(name='output_data', mode='cell',
                                 data=piezo_strain, dofs=None)
    out['total_stress'] = Struct(name='output_data', mode='cell',
                                 data=stress + piezo, dofs=None)

    return out

filename_mesh = data_dir + '/meshes/2d/special/circle_in_square.mesh'
## filename_mesh = data_dir + '/meshes/2d/special/circle_in_square_small.mesh'
## filename_mesh = data_dir + '/meshes/3d/special/cube_sphere.mesh'
## filename_mesh = data_dir + '/meshes/2d/special/cube_cylinder.mesh'

omega = 1
omega_squared = omega**2

conf_dir = os.path.dirname(__file__)
io = MeshIO.any_from_filename(filename_mesh, prefix_dir=conf_dir)
bbox, dim = io.read_bounding_box(ret_dim=True)

geom = {3 : '3_4', 2 : '2_3'}[dim]

x_left, x_right = bbox[:,0]

options = {
    'post_process_hook' : 'post_process',
}

regions = {
    'Y' : 'all',
    'Y1' : 'cells of group 1',
    'Y2' : 'cells of group 2',
    'Y2_Surface': ('r.Y1 *v r.Y2', 'facet'),
    'Left' : ('vertices in (x < %f)' % (x_left + 1e-3), 'facet'),
    'Right' : ('vertices in (x > %f)' % (x_right - 1e-3), 'facet'),
}

fields = {
    'displacement' : ('real', dim, 'Y', 1),
    'potential' : ('real', 1, 'Y', 1),
}

variables = {
    'u' : ('unknown field', 'displacement', 0),
    'v' : ('test field', 'displacement', 'u'),
    'phi' : ('unknown field', 'potential', 1),
    'psi' : ('test field', 'potential', 'phi'),
}

ebcs = {
    'u1' : ('Left', {'u.all' : 0.0}),
    'u2' : ('Right', {'u.0' : 0.1}),
    'phi' : ('Y2_Surface', {'phi.all' : 0.0}),
}

def get_inclusion_pars(ts, coor, mode=None, **kwargs):
    """TODO: implement proper 3D -> 2D transformation of constitutive
    matrices."""
    if mode == 'qp':
        _, dim = coor.shape
        sym = (dim + 1) * dim // 2

        dielectric = nm.eye(dim, dtype=nm.float64)
        # !!!
        coupling = nm.ones((dim, sym), dtype=nm.float64)
        #    coupling[0,1] = 0.2

        out = {
            # Lame coefficients in 1e+10 Pa.
            'D' : stiffness_from_lame(dim=2, lam=0.1798, mu=0.148),
            # dielectric tensor
            'dielectric' : dielectric,
            # piezoelectric coupling
            'coupling' : coupling,
            'density' : nm.array([[0.1142]]), # in 1e4 kg/m3
        }

        for key, val in out.items():
            out[key] = val[None, ...]

        return out

materials = {
    'inclusion' : (None, 'get_inclusion_pars')
}

functions = {
    'get_inclusion_pars' : (get_inclusion_pars,),
}

integrals = {
    'i' : 2,
}

equations = {
    '1' : """- %f * dw_dot.i.Y(inclusion.density, v, u)
             + dw_lin_elastic.i.Y(inclusion.D, v, u)
             - dw_piezo_coupling.i.Y2(inclusion.coupling, v, phi)
           = 0""" % omega_squared,
    '2' : """dw_piezo_coupling.i.Y2(inclusion.coupling, u, psi)
           + dw_diffusion.i.Y(inclusion.dielectric, psi, phi)
           = 0""",
}

solvers = {
    'ls' : ('ls.scipy_direct', {}),
    'newton' : ('nls.newton',
                {'i_max'      : 1,
                 'eps_a'      : 1e-10,
    }),
}
