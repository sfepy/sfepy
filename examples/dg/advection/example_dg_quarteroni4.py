"""
Based on

Antonietti, P., & Quarteroni, A. (2013). Numerical performance of discontinuous and stabilized
    continuous Galerkin methods for convection–diffusion problems.
Numerical Methods for Hyperbolic Equations, 75–85. https://doi.org/10.1201/b14172-9
"""

from examples.dg.example_dg_common import *

def define(filename_mesh=None,
           approx_order=3,

           flux=0.0,
           limit=False,

           Cw=100,
           diffusion_coef=1e-9,
           diff_scheme_name="symmetric",

           CFL=None,
           dt=None,
           ):

    example_name = "quart4"
    dim = 2

    if filename_mesh is None:
        filename_mesh = get_gen_block_mesh_hook((1., 1.), (20, 20), (.5, .5))

    velo = [3., 0.]

    angle = nm.pi / 3
    rotm = nm.array([[nm.cos(angle), nm.sin(angle)],
                     [-nm.sin(angle), nm.cos(angle)]])
    velo = nm.sum(rotm.T * nm.array(velo), axis=-1)[:, None]


    regions = {
        'Omega'     : 'all',
        'bot_left' : ('vertices in (x < 0.01) & (y < 0.5)', 'edge'),
        'top_left' : ('vertices in (x < 0.01) & (y > 0.5)', 'edge'),
        'right': ('vertices in x == 1', 'edge'),
        'top' : ('vertices in y == 1', 'edge'),
        'bottom': ('vertices in y == 0', 'edge')
    }

    fields = {
        'density': ('real', 'scalar', 'Omega', str(approx_order) + 'd', 'DG', 'legendre')
    }

    variables = {
        'u': ('unknown field', 'density', 0),
        'v': ('test field', 'density', 'u'),
    }

    integrals = {
        'i': 2 * approx_order,
    }

    dgebcs = {
        'u_bot_left' : ('bot_left', {'u.all': 1, 'grad.u.all' : (0, 0)}),
        'u_top_left': ('top_left', {'u.all': 0, 'grad.u.all': (0, 0)}),
        'u_top'  : ('top', {'u.all': 0, 'grad.u.all': (0, 0)}),
        'u_bot'  : ('bottom', {'u.all': 1, 'grad.u.all': (0, 0)}),
        'u_right': ('right', {'u.all': 0, 'grad.u.all': (0, 0)}),
    }

    materials = {
        'a'     : ({'val': [velo], '.flux': flux},),
        'D'     : ({'val': [diffusion_coef], '.Cw': Cw},),
        # 'g'     : 'source_fun'
    }

    equations = {
        'balance': """
                    + dw_s_dot_mgrad_s.i.Omega(a.val, u, v)
                    - dw_dg_advect_laxfrie_flux.i.Omega(a.flux, a.val, v, u)
                   """
                   +
                   " - dw_laplace.i.Omega(D.val, v, u) " +
                   " + dw_dg_diffusion_flux.i.Omega(D.val, u, v) " +
                   " + dw_dg_diffusion_flux.i.Omega(D.val, v, u)" +
                   " - " + str(diffusion_coef) + "* dw_dg_interior_penal.i.Omega(D.Cw, v, u)" +
                   # " + dw_volume_lvf.i.Omega(g.val, v)
                   " = 0"

    }

    solver_0 = {
        'name' : 'ls',
        'kind' : 'ls.scipy_direct',
    }

    solver_1 = {
        'name' : 'newton',
        'kind' : 'nls.newton',

        'i_max'      : 5,
        'eps_a'      : 1e-8,
        'eps_r'      : 1.0,
        'macheps'   : 1e-16,
        'lin_red'    : 1e-2,  # Linear system error < (eps_a * lin_red).
        'ls_red'     : 0.1,
        'ls_red_warp' : 0.001,
        'ls_on'      : 0.99999,
        'ls_min'     : 1e-5,
        'check'     : 0,
        'delta'     : 1e-6,
    }

    options = {
        'nls'             : 'newton',
        'ls'              : 'ls',
        'output_format'   : 'msh',
        # 'pre_process_hook': get_cfl_setup(CFL)
    }
    return locals()

globals().update(define())