import numpy as nm

from sfepy.discrete.fem import Mesh
from sfepy.discrete.fem.meshio import UserMeshIO

from sfepy.terms import register_term
from sfepy.solvers import  register_solver


# import various ICs
from sfepy.discrete.dg.my_utils.inits_consts import ghump, left_par_q, left_cos, superic, three_step_u, \
    sawtooth_q, const_q, quadr_cub

from sfepy.discrete.dg.dg_terms import AdvectDGFluxTerm


# import TSSs
from sfepy.discrete.dg.dg_tssolver import TVDRK3StepSolver, RK4StepSolver, EulerStepSolver
from sfepy.discrete.dg.dg_limiters import IdentityLimiter, Moment1DLimiter

register_term(AdvectDGFluxTerm)
register_solver(TVDRK3StepSolver)
register_solver(RK4StepSolver)
register_solver(EulerStepSolver)

example_name = "adv_1D"
dim = int(example_name[example_name.index("D") - 1])

# Setup Mesh
def get_mesh_hook(XS, XE, n_nod):


    def mesh_hook(mesh, mode):
        """
        Generate the 1D mesh.
        """
        if mode == 'read':

            coors = nm.linspace(XS, XE, n_nod).reshape((n_nod, 1))
            conn = nm.arange(n_nod, dtype=nm.int32).repeat(2)[1:-1].reshape((-1, 2))
            mat_ids = nm.zeros(n_nod - 1, dtype=nm.int32)
            descs = ['1_2']

            mesh = Mesh.from_data('laplace_1d', coors, None,
                                  [conn], [mat_ids], descs)
            return mesh

        elif mode == 'write':
            pass

    return mesh_hook

filename_mesh = UserMeshIO(get_mesh_hook(0, 1, 100))

approx_order = 1
t0 = 0.
t1 = 1.

materials = {
    'a' : ({'val': 1.0, '.flux': 0.0},),
}

regions = {
    'Omega' : 'all',
    'Gamma_Left' : ('vertices in (x < 0.00001)', 'vertex'),
    'Gamma_Right' : ('vertices in (x > 0.99999)', 'vertex'),
}

fields = {
    'density' : ('real', 'scalar', 'Omega', '1d', 'DG', 'legendre') #
}

variables = {
    'u' : ('unknown field', 'density', 0, 1),
    'v' : ('test field',    'density', 'u'),
}

ebcs = {
    'u_left' : ('Gamma_Left', {'u.all' : 0.3}),
    # 't2' : ('Gamma_Right', {'t.0' : -0.3}),
}

def get_ic(x, ic=None):
    return ghump(x - .3)

functions = {
    'get_ic' : (get_ic,)
}

ics = {
    'ic' : ('Omega', {'u.0' : 'get_ic'}),
}

integrals = {
    'i' : 2 * approx_order,
}

equations = {
    'Advection' : """
                   dw_volume_dot.i.Omega(v, u)
                   + dw_s_dot_mgrad_s.i.Omega(a.val, v, u[-1])
                   - dw_dg_advect_laxfrie_flux.i.Omega(a.val, v, u[-1]) = 0
                  """
}

solvers = {
    "tss" : ('ts.euler',
                         {"t0": t0,
                          "t1": t1,
                          "n_step": n_step, # tODO satisfy CFL condition, move to run examples script
                          'limiter' : IdentityLimiter}),
    'nls' : ('nls.newton',{} ),
    'ls'  : ('ls.scipy_direct', {})
}

options = {
    'ts' : 'tss',
    'nls' : 'newton',
    'ls' : 'ls',
    'save_times' : 100
}
