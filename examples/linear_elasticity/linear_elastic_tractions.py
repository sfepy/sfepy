r"""
Linear elasticity with pressure traction load on a surface and constrained to
one-dimensional motion.

Find :math:`\ul{u}` such that:

.. math::
    \int_{\Omega} D_{ijkl}\ e_{ij}(\ul{v}) e_{kl}(\ul{u})
    = - \int_{\Gamma_{right}} \ul{v} \cdot \ull{\sigma} \cdot \ul{n}
    \;, \quad \forall \ul{v} \;,

where

.. math::
    D_{ijkl} = \mu (\delta_{ik} \delta_{jl}+\delta_{il} \delta_{jk}) +
    \lambda \ \delta_{ij} \delta_{kl}
    \;.

and :math:`\ull{\sigma} \cdot \ul{n} = \bar{p} \ull{I} \cdot \ul{n}`
with given traction pressure :math:`\bar{p}`.
"""
from __future__ import absolute_import
import numpy as nm
from sfepy.mechanics.matcoefs import stiffness_from_lame

def linear_tension(ts, coor, mode=None, **kwargs):
    if mode == 'qp':
        val = nm.tile(1.0, (coor.shape[0], 1, 1))

        return {'val' : val}

def define():
    """Define the problem to solve."""
    from sfepy import data_dir

    filename_mesh = data_dir + '/meshes/3d/block.mesh'

    options = {
        'nls' : 'newton',
        'ls' : 'ls',
    }

    functions = {
        'linear_tension' : (linear_tension,),
    }

    fields = {
        'displacement': ('real', 3, 'Omega', 1),
    }

    materials = {
        'solid' : ({'D': stiffness_from_lame(3, lam=5.769, mu=3.846)},),
        'load' : (None, 'linear_tension')
    }

    variables = {
        'u' : ('unknown field', 'displacement', 0),
        'v' : ('test field', 'displacement', 'u'),
    }

    regions = {
        'Omega' : 'all',
        'Left' : ('vertices in (x < -4.99)', 'facet'),
        'Right' : ('vertices in (x > 4.99)', 'facet'),
    }

    ebcs = {
        'fixb' : ('Left', {'u.all' : 0.0}),
        'fixt' : ('Right', {'u.[1,2]' : 0.0}),
    }

    ##
    # Balance of forces.
    equations = {
        'elasticity' :
        """dw_lin_elastic.2.Omega( solid.D, v, u )
         = - dw_surface_ltr.2.Right( load.val, v )""",
    }

    ##
    # Solvers etc.
    solvers = {
        'ls' : ('ls.scipy_direct', {}),
        'newton' : ('nls.newton',
                    { 'i_max'      : 1,
                      'eps_a'      : 1e-10,
                      'eps_r'      : 1.0,
                      'macheps'   : 1e-16,
                      # Linear system error < (eps_a * lin_red).
                      'lin_red'    : 1e-2,
                      'ls_red'     : 0.1,
                      'ls_red_warp' : 0.001,
                      'ls_on'      : 1.1,
                      'ls_min'     : 1e-5,
                      'check'     : 0,
                      'delta'     : 1e-6,
                      })
    }

    return locals()
