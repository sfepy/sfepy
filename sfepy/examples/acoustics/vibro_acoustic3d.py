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
from sfepy import data_dir, base_dir
filename_mesh = data_dir + '/meshes/3d/acoustic_wg.vtk'

sound_speed = 343.0
wave_num = 5.5
p_inc = 300

c = sound_speed
c2 = c**2
w = wave_num * c
w2 = w**2
wc = w * c
wc2 = w * c2

regions = {
    'Omega1': 'cells of group 1',
    'Omega2': 'cells of group 2',
    'GammaIn': ('vertices of group 1', 'face'),
    'GammaOut': ('vertices of group 2', 'face'),
    'Gamma_aux': ('r.Omega1 *v r.Omega2', 'face'),
    'Gamma0_1': ('copy r.Gamma_aux', 'face', 'Omega1'),
    'Gamma0_2': ('copy r.Gamma_aux', 'face', 'Omega2'),
    'aux_Left': ('vertices in (x <  0.001)', 'face'),
    'aux_Right': ('vertices in (x > 0.299)', 'face'),
    'Gamma0_1_Left': ('r.Gamma0_1 *v r.aux_Left', 'edge'),
    'Gamma0_1_Right': ('r.Gamma0_1 *v r.aux_Right', 'edge'),
    }

fields = {
    'pressure1': ('complex', 'scalar', 'Omega1', 1),
    'pressure2': ('complex', 'scalar', 'Omega2', 1),
    'tvelocity': ('complex', 'scalar', 'Gamma0_1', 1),
    'deflection': ('complex', 'scalar', 'Gamma0_1', 1),
    }

variables = {
    'p1': ('unknown field', 'pressure1', 0),
    'q1': ('test field', 'pressure1', 'p1'),
    'p2': ('unknown field', 'pressure2', 1),
    'q2': ('test field', 'pressure2', 'p2'),
    'g0': ('unknown field', 'tvelocity', 2),
    'f0': ('test field', 'tvelocity', 'g0'),
    'w': ('unknown field', 'deflection', 3),
    'z': ('test field', 'deflection', 'w'),
    }

ebcs = {
    'fixed_l': ('Gamma0_1_Left', {'w.0': 0.0}),
    'fixed_r': ('Gamma0_1_Right', {'w.0': 0.0}),
    }

options = {
    'split_results_by': 'region',
    'output_dir': 'output',
    }

functions = {
    }

materials = {
    'ac' : ({'F': -2.064e+00, 'c': -1.064e+00}, ),
    }

equations = {
    'eq_1' : """
        %e * dw_laplace.5.Omega1(q1, p1)
      + %e * dw_laplace.5.Omega2(q2, p2)
      - %e * dw_dot.5.Omega1(q1, p1)
      - %e * dw_dot.5.Omega2(q2, p2)
      + %s * dw_dot.5.GammaIn(q1, p1)
      + %s * dw_dot.5.GammaOut(q2, p2)
      - %s * dw_dot.5.Gamma0_1(q1, g0)
      + %s * dw_dot.5.Gamma0_2(q2, tr(g0))
      = %s * dw_integrate.5.GammaIn(q1)"""\
        % (c2, c2, w2, w2,
           1j * wc, 1j * wc,
           1j * wc2, 1j * wc2,
           2j * wc * p_inc),
    'eq_2' : """
      - %s * dw_dot.5.Gamma0_1(f0, p1)
      + %s * dw_dot.5.Gamma0_1(f0, tr(p2))
      - %e * dw_dot.5.Gamma0_1(ac.F, f0, g0)
      + %e * dw_dot.5.Gamma0_1(ac.c, f0, w)
      = 0"""\
        % (1j * w, 1j * w, w2, w2),
    }

solvers = {
    'ls': ('ls.cm_pb',
           {'others': [base_dir
                       + '/examples/acoustics/vibro_acoustic3d_mid.py'],
            'coupling_variables': ['g0', 'w'],
            }),
    'nls': ('nls.newton', {
        'i_max' : 1,
        'eps_a' : 1e-6,
        'eps_r' : 1e-6,
    })
}
