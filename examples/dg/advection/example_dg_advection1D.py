from examples.dg.example_dg_common import *

example_name = "adv_1D"
dim = int(example_name[example_name.index("D") - 1])


# Setup Mesh
def get_1Dmesh_hook(XS, XE, n_nod):
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


filename_mesh = UserMeshIO(get_1Dmesh_hook(0, 1, 100))

approx_order = 2
t0 = 0.
t1 = 1.
CFL = .5

materials = {
    'a': ({'val': 1.0, '.flux': 0.0},),
}

regions = {
    'Omega'      : 'all',
    'Gamma_Left' : ('vertices in (x < 0.0001)', 'vertex'),
    'Gamma_Right': ('vertices in (x > 0.9995)', 'vertex'),
}

fields = {
    'density': ('real', 'scalar', 'Omega', str(approx_order) + 'd', 'DG', 'legendre')
}

variables = {
    'u': ('unknown field', 'density', 0, 1),
    'v': ('test field', 'density', 'u'),
}


def left_sin(ts, coor, bc, problem, **kwargs):
    return nm.sin(2 * nm.pi * ts.time)


def get_ic(x, ic=None):
    return ghump(x - .6)


from sfepy.discrete.fem.periodic import match_y_line

functions = {
    'get_ic'      : (get_ic,),
    'bc_fun'      : (left_sin,),
    'match_y_line': (match_y_line,)
}

# ebcs = {
#     'u_left' : ('Gamma_Left', {'u.all' : .5}),
#     # 'u_righ' : ('Gamma_Right', {'u.all' : -0.3}),
# }

epbc_1 = {
    'name'  : 'u_rl',
    'region': ['Gamma_Left', 'Gamma_Right'],
    'dofs'  : {'u.all': 'u.all'},
    'match' : 'match_y_line',
}

ics = {
    'ic': ('Omega', {'u.0': 'get_ic'}),
}

integrals = {
    'i': 2 * approx_order,
}

equations = {
    'Advection': """
                   dw_volume_dot.i.Omega(v, u)
                   + dw_s_dot_mgrad_s.i.Omega(a.val, u[-1], v)
                   - dw_dg_advect_laxfrie_flux.i.Omega(a.val, v, u[-1]) = 0
                  """
}

solvers = {
    "tss": ('ts.euler',
            {"t0"     : t0,
             "t1"     : t1,
             'limiter': IdentityLimiter,
             'verbose': True}),
    'nls': ('nls.newton', {}),
    'ls' : ('ls.scipy_direct', {})
}

options = {
    'ts'              : 'tss',
    'nls'             : 'newton',
    'ls'              : 'ls',
    'save_times'      : 100,
    'active_only'     : False,
    'pre_process_hook': get_cfl_setup(CFL),
    'output_format'   : "vtk"
}
