r"""
Thermo-elasticity with a given temperature distribution.

Uses `dw_biot` term with an isotropic coefficient for thermo-elastic coupling.

For given body temperature :math:`T` and background temperature
:math:`T_0` find :math:`\ul{u}` such that:

.. math::
    \int_{\Omega} D_{ijkl}\ e_{ij}(\ul{v}) e_{kl}(\ul{u})
    - \int_{\Omega}  (T - T_0)\ \alpha_{ij} e_{ij}(\ul{v})
    = 0
    \;, \quad \forall \ul{v} \;,

where

.. math::
    D_{ijkl} = \mu (\delta_{ik} \delta_{jl}+\delta_{il} \delta_{jk}) +
    \lambda \ \delta_{ij} \delta_{kl}
    \;, \\

    \alpha_{ij} = (3 \lambda + 2 \mu) \alpha \delta_{ij}

and :math:`\alpha` is the thermal expansion coefficient.
"""
from __future__ import absolute_import
import numpy as np

from sfepy.base.base import Struct
from sfepy.mechanics.matcoefs import stiffness_from_lame
from sfepy.mechanics.tensors import get_von_mises_stress
from sfepy import data_dir

# Material parameters.
lam = 10.0
mu = 5.0
thermal_expandability = 1.25e-5
T0 = 20.0 # Background temperature.

filename_mesh = data_dir + '/meshes/3d/block.mesh'

def get_temperature_load(ts, coors, region=None, **kwargs):
    """
    Temperature load depends on the `x` coordinate.
    """
    x = coors[:, 0]
    return (x - x.min())**2 - T0

def post_process(out, pb, state, extend=False):
    """
    Compute derived quantities: strain, stresses. Store also the loading
    temperature.
    """
    ev = pb.evaluate

    strain = ev('ev_cauchy_strain.2.Omega( u )', mode='el_avg')
    out['cauchy_strain'] = Struct(name='output_data',
                                  mode='cell', data=strain,
                                  dofs=None)

    e_stress = ev('ev_cauchy_stress.2.Omega( solid.D, u )', mode='el_avg')
    out['elastic_stress'] = Struct(name='output_data',
                                   mode='cell', data=e_stress,
                                   dofs=None)

    t_stress = ev('ev_biot_stress.2.Omega( solid.alpha, T )', mode='el_avg')
    out['thermal_stress'] = Struct(name='output_data',
                                   mode='cell', data=t_stress,
                                   dofs=None)

    out['total_stress'] = Struct(name='output_data',
                                 mode='cell', data=e_stress + t_stress,
                                 dofs=None)

    out['von_mises_stress'] = aux = out['total_stress'].copy()
    vms = get_von_mises_stress(aux.data.squeeze())
    vms.shape = (vms.shape[0], 1, 1, 1)
    out['von_mises_stress'].data = vms

    val = pb.get_variables()['T']()
    val.shape = (val.shape[0], 1)
    out['T'] = Struct(name='output_data',
                      mode='vertex', data=val + T0,
                      dofs=None)
    return out

options = {
    'post_process_hook' : 'post_process',

    'nls' : 'newton',
    'ls' : 'ls',
}

functions = {
    'get_temperature_load' : (get_temperature_load,),
}

regions = {
    'Omega' : 'all',
    'Left' : ('vertices in (x < -4.99)', 'facet'),
}

fields = {
    'displacement': ('real', 3, 'Omega', 1),
    'temperature': ('real', 1, 'Omega', 1),
}

variables = {
    'u' : ('unknown field', 'displacement', 0),
    'v' : ('test field', 'displacement', 'u'),
    'T' : ('parameter field', 'temperature',
           {'setter' : 'get_temperature_load'}),
}

ebcs = {
    'fix_u' : ('Left', {'u.all' : 0.0}),
}

eye_sym = np.array([[1], [1], [1], [0], [0], [0]], dtype=np.float64)
materials = {
    'solid' : ({
        'D' : stiffness_from_lame(3, lam=lam, mu=mu),
        'alpha' : (3.0 * lam + 2.0 * mu) * thermal_expandability * eye_sym
    },),
}

equations = {
    'balance_of_forces' :
    """dw_lin_elastic.2.Omega( solid.D, v, u )
     - dw_biot.2.Omega( solid.alpha, v, T )
     = 0""",
}

solvers = {
    'ls' : ('ls.scipy_direct', {}),
    'newton' : ('nls.newton', {
        'i_max'      : 1,
        'eps_a'      : 1e-10,
    }),
}
