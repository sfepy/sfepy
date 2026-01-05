import numpy as np

# Material
E = 70e9
nu = 0.3

# Geometry
r_inner, r_outer = 0.1, 0.2
height = 0.3

# Manufactured solution
a, b = 1e-3, 2e-3  # u_r = a*r, u_z = b*z

def define():
    import tempfile
    from sfepy.discrete.fem import Mesh
    
    # Lame parameters
    lam = E * nu / ((1 + nu) * (1 - 2*nu))
    mu = E / (2 * (1 + nu))
    
    # Create mesh
    nr, nz = 9, 11
    r = np.linspace(r_inner, r_outer, nr)
    z = np.linspace(0, height, nz)
    
    nodes = []
    for iz in range(nz):
        for ir in range(nr):
            nodes.append([r[ir], z[iz]])
    coors = np.array(nodes)
    
    conn = []
    for iz in range(nz - 1):
        for ir in range(nr - 1):
            n0 = iz * nr + ir
            conn.append([n0, n0 + 1, n0 + nr + 1, n0 + nr])
    conn = np.array(conn)
    
    mesh = Mesh.from_data('annulus', coors, None, [conn],
                         [np.zeros(len(conn), dtype=np.int32)], ['2_4'])
    
    filename_mesh = tempfile.mktemp(suffix='.vtk', prefix='axisym_')
    mesh.write(filename_mesh, io='auto')
    
    # Material coefficient functions
    # The key idea: encode 2πr, 1/r, and material properties in coefficients
    
    def material_rr(ts, coors, mode=None, **kwargs):
        """Coefficient for ∂u_r/∂r term: r * [[(λ+2μ), 0], [0, μ]]"""
        if mode == 'qp':
            r = coors[:, 0]
            K = np.zeros((coors.shape[0], 2, 2), dtype=np.float64)
            K[:, 0, 0] = (lam + 2*mu) * r
            K[:, 1, 1] = mu * r
            return {'K': K}
        return None
    
    def material_zz(ts, coors, mode=None, **kwargs):
        """Coefficient for ∂u_z/∂z term: r * [[μ, 0], [0, (λ+2μ)]]"""
        if mode == 'qp':
            r = coors[:, 0]
            K = np.zeros((coors.shape[0], 2, 2), dtype=np.float64)
            K[:, 0, 0] = mu * r
            K[:, 1, 1] = (lam + 2*mu) * r
            return {'K': K}
        return None
    
    def material_poisson_rr(ts, coors, mode=None, **kwargs):
        """Poisson coupling for u_r equation: r * [[0, μ], [λ, 0]]"""
        if mode == 'qp':
            r = coors[:, 0]
            K = np.zeros((coors.shape[0], 2, 2), dtype=np.float64)
            K[:, 0, 1] = mu * r
            K[:, 1, 0] = lam * r
            return {'K': K}
        return None
    
    def material_poisson_zz(ts, coors, mode=None, **kwargs):
        """Poisson coupling for u_z equation: r * [[0, λ], [μ, 0]]"""
        if mode == 'qp':
            r = coors[:, 0]
            K = np.zeros((coors.shape[0], 2, 2), dtype=np.float64)
            K[:, 0, 1] = lam * r
            K[:, 1, 0] = mu * r
            return {'K': K}
        return None
    
    def material_hoop_idx(ts, coors, mode=None, **kwargs):
        """Hoop strain coupling coefficient: scalar λ"""
        if mode == 'qp':
            val = lam * np.ones((coors.shape[0], 1, 1), dtype=np.float64)
            return {'val': val}
        return None
    
    def material_hoop0(ts, coors, mode=None, **kwargs):
        """Hoop strain volume term: c = (λ+2μ)/r"""
        if mode == 'qp':
            r = coors[:, 0]
            r_safe = np.maximum(r, 1e-10)
            c = (lam + 2*mu) / r_safe
            return {'c': c.reshape(-1, 1, 1)}
        return None
    
    regions = {
        'Omega': 'all',
        'Left': ('vertices in (x < %.10f)' % (r_inner + 1e-6), 'facet'),
        'Right': ('vertices in (x > %.10f)' % (r_outer - 1e-6), 'facet'),
        'Bottom': ('vertices in (y < %.10f)' % 1e-6, 'facet'),
        'Top': ('vertices in (y > %.10f)' % (height - 1e-6), 'facet'),
    }
    
    materials = {
        'm_rr': 'material_rr',
        'm_zz': 'material_zz',
        'm_prr': 'material_poisson_rr',
        'm_pzz': 'material_poisson_zz',
        'm_hoop_idx': 'material_hoop_idx',
        'm_hoop0': 'material_hoop0',
    }
    
    fields = {
    'ur': ('real', 1, 'Omega', 1),
    'uz': ('real', 1, 'Omega', 1),
}
    
    integrals = {'i': 2}
    
    variables = {
    'ur': ('unknown field', 'ur', 0),
    'vr': ('test field',    'ur', 'ur'),
    'uz': ('unknown field', 'uz', 1),
    'vz': ('test field',    'uz', 'uz'),
}
    
    functions = {
        'material_rr': (material_rr,),
        'material_zz': (material_zz,),
        'material_poisson_rr': (material_poisson_rr,),
        'material_poisson_zz': (material_poisson_zz,),
        'material_hoop_idx': (material_hoop_idx,),
        'material_hoop0': (material_hoop0,),
        'bc_ur_left': (lambda ts, coors, **kw: a * r_inner * np.ones(len(coors)),),
        'bc_ur_right': (lambda ts, coors, **kw: a * r_outer * np.ones(len(coors)),),
        'bc_ur_bottom': (lambda ts, coors, **kw: a * coors[:, 0],),
        'bc_ur_top': (lambda ts, coors, **kw: a * coors[:, 0],),
        'bc_uz_left': (lambda ts, coors, **kw: b * coors[:, 1],),
        'bc_uz_right': (lambda ts, coors, **kw: b * coors[:, 1],),
        'bc_uz_bottom': (lambda ts, coors, **kw: np.zeros(len(coors)),),
        'bc_uz_top': (lambda ts, coors, **kw: b * height * np.ones(len(coors)),),
    }
    
    # Boundary conditions
    ebcs = {
        'fix_ur_left': ('Left', {'u.0': 'bc_ur_left'}),
        'fix_ur_right': ('Right', {'u.0': 'bc_ur_right'}),
        'fix_ur_bottom': ('Bottom', {'u.0': 'bc_ur_bottom'}),
        'fix_ur_top': ('Top', {'u.0': 'bc_ur_top'}),
        'fix_uz_left': ('Left', {'u.1': 'bc_uz_left'}),
        'fix_uz_right': ('Right', {'u.1': 'bc_uz_right'}),
        'fix_uz_bottom': ('Bottom', {'u.1': 'bc_uz_bottom'}),
        'fix_uz_top': ('Top', {'u.1': 'bc_uz_top'}),
    }
    
    # Equations decomposing the axisymmetric weak form
    equations = {
        'balance_r': r"""
            dw_diffusion.i.Omega(m_rr.K, vr, ur)
        + dw_diffusion.i.Omega(m_prr.K, vr, uz)
        - dw_s_dot_grad_i_s.i.Omega(m_hoop_idx.val, vr, uz)
        + dw_volume_dot.i.Omega(m_hoop0.c, vr, ur)
        = 0
        """,
        'balance_z': r"""
            dw_diffusion.i.Omega(m_zz.K, vz, uz)
        + dw_diffusion.i.Omega(m_pzz.K, vz, ur)
        + dw_s_dot_grad_i_s.i.Omega(m_hoop_idx.val, vz, ur)
        = 0
        """,
    }

    
    solvers = {
        'ls': ('ls.scipy_direct', {}),
        'newton': ('nls.newton', {'i_max': 1, 'eps_a': 1e-10}),
    }
    
    options = {
        'nls': 'newton',
        'ls': 'ls',
        'post_process_hook': 'post_process',
    }
    
    return locals()

def post_process(out, problem, variables, extend=False):
    """Calculate error against manufactured solution."""
    ur = variables['ur']
    uz = variables['uz']
    
    coors = problem.domain.get_mesh_coors()
    r, z = coors[:, 0], coors[:, 1]
    
    # Exact solution
    u_r_exact = a * r
    u_z_exact = b * z
    
    # Computed
    u_r_val = ur().ravel()
    u_z_val = uz().ravel()
    
    # Errors
    err_r = np.abs(u_r_val - u_r_exact)
    err_z = np.abs(u_z_val - u_z_exact)
    err_total = np.sqrt(err_r**2 + err_z**2)
    
    max_err = np.max(err_total)
    
    print()
    print('='*70)
    print('AXISYMMETRIC ELASTICITY VIA TERM COMPOSITION')
    print('='*70)
    print('Approach: Decompose weak form into existing SfePy terms')
    print('  - dw_laplace: for gradient terms (∂u/∂r, ∂u/∂z)')
    print('  - dw_volume_dot: for hoop strain (u_r/r)')
    print('  - Spatial coefficients: 2πr weight + material properties')
    print()
    print('Manufactured solution: u_r = {:.1e}·r, u_z = {:.1e}·z'.format(a, b))
    print('-'*70)
    print('Maximum displacement error: {:.3e}'.format(max_err))
    print()
    if max_err < 1e-8:
        print('✓ PASS: Error < 1e-8 (expected for linear elements)')
    elif max_err < 1e-4:
        print('⚠ PARTIAL: Acceptable but not optimal')
    else:
        print('✗ FAIL: Error too large - check implementation')
    print('='*70)
    print()
    
    # Note: We don't store outputs because we've already printed the results
    # and storing raw numpy arrays causes issues with SfePy's output system
    
    return out
