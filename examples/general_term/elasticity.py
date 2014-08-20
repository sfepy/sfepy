r"""
This example shows a functionality of general term for elasticity with a very
simple mesh.

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
"""
from sfepy.mechanics.matcoefs import lame_from_youngpoisson
from sfepy.discrete.fem.utils import refine_mesh
from sfepy import data_dir

filename_mesh = data_dir + '/meshes/2d/simple_square.mesh'

young = 2000.0 # Young's modulus [MPa]
poisson = 0.4  # Poisson's ratio
from sfepy.mechanics.matcoefs import ElasticConstants, ElasticTensor
elas = ElasticConstants(young=2000., poisson=.4)
bulk, mu = elas.get(['bulk', 'mu'])
elasTensor = ElasticTensor(bulk=bulk, mu=mu, plane='strain')

materials = {
    'Asphalt' : ({
        'lam' : lame_from_youngpoisson(young, poisson)[0],
        'mu' : lame_from_youngpoisson(young, poisson)[1],
        'Voight' : elasTensor.voight,
        'Mandel' : elasTensor.mandel,
    },),
    'Load' : ({'.val' : [0.0, -1000.0]},),
}

regions = {
    'Omega' : 'all',
    'Omega1' : 'cells of group 1',
    'Omega2' : 'cells of group 2',
    'Left' : ('vertices in (x < 0.0001)', 'facet'),
    'Right' : ('vertices in (x > 0.9999)', 'facet'),
    'Bottom' : ('vertices in (y < 0.0001)', 'facet'),
    'Top' : ('vertices in (y > 0.9999)', 'facet'),
    'Vertex' : ('vertex 3', 'vertex'),
}

fields = {
    'displacement': ('real', 'vector', 'Omega', 1),
}

variables = {
    'u' : ('unknown field', 'displacement', 0),
    'v' : ('test field', 'displacement', 'u'),
}

equations = {
   'balance_of_forces' :
   """intFE.2.Omega.einsum('ij,j,i', Asphalt.Mandel, v.grad.sym.Man,
                           u.grad.sym.Man)
      = dw_point_load.0.Vertex(Load.val, v)""",
}


ebcs = {
    'XSym' : ('Bottom', {'u.1' : 0.0}),
    'YSym' : ('Left', {'u.0' : 0.0}),
}

solvers = {
    'ls' : ('ls.scipy_direct', {}),
    'newton' : ('nls.newton', {
        'i_max' : 1,
        'eps_a' : 1e-6,
        'problem' : 'nonlinear'
    }),
}

options = {'post_process_hook_final': 'post_process_hook_final'}

def post_process_hook_final(pb, state):
    ev = pb.evaluate
    energy_norm = []
    energy_norm.append(ev("dw_lin_elastic.2.Omega(Asphalt.Voight, u, u)",
                          mode='eval'))
    eval_str = """intFE.2.Omega.einsum('ij,j,i', Asphalt.Mandel,
                                        u.grad.sym.Man, u.grad.sym.Man)"""
    energy_norm.append(ev(eval_str, mode='eval'))
    print "Energetic norms for dw_lin_elastic term (%s) and intFE term (%s)." \
        % (str(energy_norm[0]), str(energy_norm[1]))

