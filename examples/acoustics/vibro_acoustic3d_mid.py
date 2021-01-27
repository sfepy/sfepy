r"""
Vibro-acoustic problem

3D acoustic domain with 2D perforated deforming interface.

*Master problem*: defined in 3D acoustic domain (``vibro_acoustic3d.py``)

*Slave subproblem*: 2D perforated interface (``vibro_acoustic3d_mid.py``)

Master 3D problem - find :math:`p` (acoustic pressure)
and :math:`g` (transversal acoustic velocity) such that:

.. math::
    c^2 \int_{\Omega} \nabla q \cdot \nabla p
    - \omega^2 \int_{\Omega} q p
    + i \omega c \int_{\Gamma_{in}} q p
    + i \omega c \int_{\Gamma_{out}} q p
    - i \omega c^2 \int_{\Gamma_0} (q^+ - q^-) g
    = 2i \omega c \int_{\Gamma_{in}} q \bar{p}
    \;, \quad \forall q \;,

    - i \omega \int_{\Gamma_0} f (p^+ - p^-)
    - \omega^2 \int_{\Gamma_0} F f g
    + \omega^2 \int_{\Gamma_0} C f w
    = 0
    \;, \quad \forall f \;,

Slave 2D subproblem - find :math:`w` (plate deflection)
and :math:`\ul{\theta}` (rotation) such that:

.. math::
    \omega^2 \int_{\Gamma_0} C z g
    - \omega^2 \int_{\Gamma_0} S z w
    + \int_{\Gamma_0} \nabla z \cdot \ull{G} \cdot \nabla w
    - \int_{\Gamma_0} \ul{\theta} \cdot \ull{G} \cdot \nabla z
    = 0
    \;, \quad \forall z \;,

    - \omega^2 \int_{\Gamma_0} R\, \ul{\nu} \cdot \ul{\theta}
    + \int_{\Gamma_0} D_{ijkl} e_{ij}(\ul{\nu}) e_{kl}(\ul{\theta})
    - \int_{\Gamma_0} \ul{\nu} \cdot \ull{G} \cdot \nabla w
    + \int_{\Gamma_0} \ul{\nu} \cdot \ull{G} \cdot \ul{\theta}
    = 0
    \;, \quad \forall \ul{\nu} \;,
"""
from __future__ import absolute_import
import numpy as nm
from sfepy.mechanics.matcoefs import stiffness_from_lame

filename_mesh = '../../meshes/2d/acoustic_wg_mid.vtk'

sound_speed = 343.0
wave_num = 5.5
thickness = 0.01

c = sound_speed
c2 = c**2
w = wave_num * c
w2 = w**2
wc = w * c
wc2 = w * c2

regions = {
    'Gamma0': 'all',
    'Left': ('vertices in (x < 0.001)', 'facet'),
    'Right': ('vertices in (x > 0.299)', 'facet'),
    }

fields = {
    'deflection': ('complex', 'scalar', 'Gamma0', 1),
    'rotation': ('complex', 'vector', 'Gamma0', 1),
    'tvelocity': ('complex', 'scalar', 'Gamma0', 1),
    }

variables = {
    'w': ('unknown field', 'deflection'),
    'z': ('test field', 'deflection', 'w'),
    'theta': ('unknown field', 'rotation'),
    'nu': ('test field', 'rotation', 'theta'),
    'g0': ('unknown field', 'tvelocity'),
    'f0': ('test field', 'tvelocity', 'g0'),
    }

ebcs = {
    'fixed_l': ('Left', {'w.0': 0.0, 'theta.all': 0.0}),
    'fixed_r': ('Right', {'w.0': 0.0, 'theta.all': 0.0}),
    }

options = {
    }

materials = {
    'ac' : ({'c': -1.064e+00, 'T': 9.202e-01,
            'hG': thickness * 4.5e10 * nm.eye(2),
            'hR': thickness * 0.71,
            'h3R': thickness**3 / 3.0 * 0.71,
            'h3C': thickness**3 / 3.0 * stiffness_from_lame(2, 1e1, 1e0)}, ),
    }

equations = {
    'eq_3': """
        %e * dw_volume_dot.5.Gamma0(ac.c, z, g0)
      - %e * dw_volume_dot.5.Gamma0(ac.T, z, w)
      - %e * dw_volume_dot.5.Gamma0(ac.hR, z, w)
           + dw_diffusion.5.Gamma0(ac.hG, z, w)
           - dw_v_dot_grad_s.5.Gamma0(ac.hG, theta, z)
      = 0"""\
        % (w2, w2, w2),
    'eq_4': """
      - %e * dw_volume_dot.5.Gamma0(ac.h3R, nu, theta)
           + dw_lin_elastic.5.Gamma0(ac.h3C, nu, theta)
           - dw_v_dot_grad_s.5.Gamma0(ac.hG, nu, w)
           + dw_volume_dot.5.Gamma0(ac.hG, nu, theta)
      = 0"""\
        % (w2, ),
    }

solvers = {
    'ls' : ('ls.scipy_direct', {}),
    'newton' : ('nls.newton', {
        'i_max' : 1,
        'eps_a' : 1e-4,
        'eps_r' : 1e-4,
    })
}
