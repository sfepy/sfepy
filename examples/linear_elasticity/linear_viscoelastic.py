r"""
Linear viscoelasticity with pressure traction load on surface
:math:`\Gamma_{right}` and constrained to one-dimensional motion.

The fading memory terms require an unloaded initial configuration, so
the load starts in the second time step.

This example uses exponential fading memory kernel
:math:`\Hcal_{ijkl}(t) = \Hcal_{ijkl}(0) e^{-d t}` with decay
:math:`d`. Two equation kinds are supported - 'th' and 'eth'. In 'th'
mode the tabulated kernel is linearly interpolated to required times
using :func:`interp_conv_mat()`. In 'eth' mode, the computation is exact
for exponential kernels.

Find :math:`\ul{u}` such that:

.. math::
    \int_{\Omega} D_{ijkl}\ e_{ij}(\ul{v}) e_{kl}(\ul{u}) \\
    + \int_{\Omega} \left [\int_0^t
    \Hcal_{ijkl}(t-\tau)\,e_{kl}(\pdiff{\ul{u}}{\tau}(\tau)) \difd{\tau}
    \right]\,e_{ij}(\ul{v}) \\
    = - \int_{\Gamma_{right}} \ul{v} \cdot \ull{\sigma} \cdot \ul{n}
    \;, \quad \forall \ul{v} \;,

where

.. math::
    D_{ijkl} = \mu (\delta_{ik} \delta_{jl}+\delta_{il} \delta_{jk}) +
    \lambda \ \delta_{ij} \delta_{kl}
    \;,

:math:`\Hcal_{ijkl}(0)` has the same structure as :math:`D_{ijkl}` and
:math:`\ull{\sigma} \cdot \ul{n} = \bar{p} \ull{I} \cdot \ul{n}` with
given traction pressure :math:`\bar{p}`.
"""
import numpy as nm

from sfepy.base.base import output
from sfepy.mechanics.matcoefs import stiffness_tensor_lame
from sfepy.homogenization.utils import interp_conv_mat
from sfepy import data_dir

def linear_tension(ts, coors, mode=None, **kwargs):
    if mode == 'qp':
        val = 1.0 * (ts.step > 0)
        output('load:', val)
        val = nm.tile(val, (coors.shape[0], 1, 1))

        return {'val' : val}

def get_exp_fading_kernel(coef0, decay, times):
    val = coef0[None, ...] * nm.exp(-decay * times[:, None, None])
    return val

def get_th_pars(ts, coors, mode=None, times=None, kernel=None, **kwargs):
    out = {}

    if mode == 'special':
        out['H'] = interp_conv_mat(kernel, ts, times)

    elif mode == 'qp':
        out['H0'] = kernel[0]
        out['Hd'] = kernel[1, 0, 0] / kernel[0, 0, 0]

        for key, val in out.iteritems():
            out[key] = nm.tile(val, (coors.shape[0], 1, 1))

    return out

filename_mesh = data_dir + '/meshes/3d/block.mesh'

# Time stepping times.
t0 = 0.0
t1 = 20.0
n_step = 21

# Fading memory times.
f_t0 = 0.0
f_t1 = 10.0
f_n_step = 11

decay = 0.5
mode = 'eth'

times = nm.linspace(f_t0, f_t1, f_n_step)
kernel = get_exp_fading_kernel(stiffness_tensor_lame(3, lam=1.0, mu=1.0),
                               decay, times)

dt = (t1 - t0) / (n_step - 1)
fading_memory_length = min(int((f_t1 - f_t0) / dt), n_step)

output('fading memory length:', fading_memory_length)

options = {
    'ts' : 'ts',
    'nls' : 'newton',
    'ls' : 'ls',
}

functions = {
    'linear_tension' : (linear_tension,),
    'get_pars' : (lambda ts, coors, mode=None, **kwargs:
                  get_th_pars(ts, coors, mode, times=times, kernel=kernel,
                              **kwargs),),
}

fields = {
    'displacement': ('real', 3, 'Omega', 1),
}

materials = {
    'solid' : ({
        'lam' : 5.769,
        'mu' : 3.846,
    },),
    'th' : 'get_pars',
    'load' : 'linear_tension',
}

variables = {
    'u' : ('unknown field', 'displacement', 0, fading_memory_length),
    'v' : ('test field', 'displacement', 'u'),
}

regions = {
    'Omega' : ('all', {}),
    'Left' : ('nodes in (x < -4.99)', {}),
    'Right' : ('nodes in (x > 4.99)', {}),
}

ebcs = {
    'fixb' : ('Left', {'u.all' : 0.0}),
    'fixt' : ('Right', {'u.[1,2]' : 0.0}),
}

if mode == 'th':
    # General form with tabulated kernel.
    equations = {
        'elasticity' :
        """dw_lin_elastic_iso.2.Omega( solid.lam, solid.mu, v, u )
         + dw_lin_elastic_th.2.Omega( ts, th.H, v, du/dt )
         = - dw_surface_ltr.2.Right( load.val, v )""",
    }

else:
    # Fast form that is exact for exponential kernels.
    equations = {
        'elasticity' :
        """dw_lin_elastic_iso.2.Omega( solid.lam, solid.mu, v, u )
         + dw_lin_elastic_eth.2.Omega( ts, th.H0, th.Hd, v, du/dt )
         = - dw_surface_ltr.2.Right( load.val, v )""",
    }

solvers = {
    'ls' : ('ls.scipy_direct', {}),
    'newton' : ('nls.newton',
                { 'i_max' : 1,
                  'eps_a' : 1e-10,
                  'problem' : 'nonlinear'}),
    'ts' : ('ts.simple', {
        't0' : t0,
        't1' : t1,
        'dt' : None,
        'n_step' : n_step,
        'quasistatic' : True,
    }),
}
