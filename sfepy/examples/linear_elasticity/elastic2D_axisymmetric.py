r"""
Axisymmetric linear elasticity using term composition approach (No custom C code).

================================================================================
                    完整弱形式展开推导（Weak Form Derivation）
================================================================================

1. 强形式（Strong Form）
------------------------
轴对称问题：在柱坐标系(r, θ, z)中，假设解与θ无关（轴对称）

位移场：u = (u_r(r,z), 0, u_z(r,z))

应变分量：
    ε_rr = ∂u_r/∂r                      （径向正应变）
    ε_zz = ∂u_z/∂z                      （轴向正应变）
    ε_θθ = u_r/r                        （环向正应变，几何项！）
    ε_rz = (∂u_r/∂z + ∂u_z/∂r)/2       （剪切应变）

应变向量：ε = [ε_rr, ε_zz, ε_θθ, ε_rz]^T ∈ ℝ⁴

本构关系（各向同性）：
    σ = D ε
    
其中刚度矩阵 D（4×4）：
    D = ⎡ λ+2μ    λ       λ      0  ⎤
        ⎢  λ     λ+2μ    λ      0  ⎥
        ⎢  λ      λ     λ+2μ    0  ⎥
        ⎣  0      0       0      μ  ⎦

平衡方程（无体力）：
    ∂σ_rr/∂r + ∂σ_rz/∂z + (σ_rr - σ_θθ)/r = 0   （径向）
    ∂σ_rz/∂r + ∂σ_zz/∂z + σ_rz/r = 0             （轴向）


2. 虚功原理（Principle of Virtual Work）
-----------------------------------------
对任意虚位移 v = (v_r, v_z)，内部虚功等于外部虚功：

    δW_int = δW_ext

对于轴对称问题，体积微元 dV = 2πr dr dz，因此：

    ∫_V σ:ε(v) dV = ∫_{Ω_rz} 2πr [σ:ε(v)] dr dz

其中 Ω_rz 是 r-z 平面上的截面区域。


3. 虚应变（Virtual Strain）
----------------------------
虚位移 v 对应的虚应变：
    ε(v) = [∂v_r/∂r, ∂v_z/∂z, v_r/r, (∂v_r/∂z + ∂v_z/∂r)/2]^T


4. 应力应变内积展开
--------------------
σ:ε(v) = σ_rr ε_rr(v) + σ_zz ε_zz(v) + σ_θθ ε_θθ(v) + 2σ_rz ε_rz(v)

代入 σ = Dε：

σ_rr = (λ+2μ)ε_rr + λε_zz + λε_θθ
σ_zz = λε_rr + (λ+2μ)ε_zz + λε_θθ
σ_θθ = λε_rr + λε_zz + (λ+2μ)ε_θθ
σ_rz = μ·2ε_rz = μ(∂u_r/∂z + ∂u_z/∂r)


5. 弱形式（逐项展开）
----------------------
∫_{Ω} 2πr [σ:ε(v)] dΩ = 0

= ∫ 2πr [(λ+2μ)ε_rr·ε_rr(v) + λε_zz·ε_zz(v) + λε_θθ·ε_θθ(v)     （正应变能）
         + λε_rr·ε_zz(v) + λε_rr·ε_θθ(v)                        （泊松耦合1）
         + λε_zz·ε_rr(v) + λε_zz·ε_θθ(v)                        （泊松耦合2）
         + λε_θθ·ε_rr(v) + λε_θθ·ε_zz(v) + (λ+2μ)ε_θθ·ε_θθ(v)  （环向耦合）
         + 2μ·2ε_rz·ε_rz(v)] dΩ                                 （剪切能）

代入应变定义：

= ∫ 2πr [(λ+2μ)·∂u_r/∂r·∂v_r/∂r                               （项1：径向梯度）
       + (λ+2μ)·∂u_z/∂z·∂v_z/∂z                               （项2：轴向梯度）
       + λ·∂u_r/∂r·∂v_z/∂z                                     （项3：泊松耦合 rr-zz）
       + λ·∂u_z/∂z·∂v_r/∂r                                     （项4：泊松耦合 zz-rr）
       + λ·∂u_r/∂r·v_r/r                                       （项5：泊松耦合 rr-θθ）
       + λ·u_r/r·∂v_r/∂r                                       （项6：泊松耦合 θθ-rr）
       + λ·∂u_z/∂z·v_r/r                                       （项7：泊松耦合 zz-θθ）
       + λ·u_r/r·∂v_z/∂z                                       （项8：泊松耦合 θθ-zz）
       + (λ+2μ)·u_r/r·v_r/r                                    （项9：环向能量）
       + 2μ(∂u_r/∂z + ∂u_z/∂r)(∂v_r/∂z + ∂v_z/∂r)] dΩ        （项10：剪切能）


6. 合并同类项
--------------
径向方程（对 v_r 的系数）：

∫ 2πr [(λ+2μ)·∂u_r/∂r·∂v_r/∂r              ← dw_laplace(2πr(λ+2μ), v_r, u_r)
      + λ·∂u_z/∂z·∂v_r/∂r                   ← dw_laplace(2πrλ, v_r, u_z)
      + λ·u_r/r·∂v_r/∂r                     ← 合并到上面
      + (λ+2μ-λ-λ)·u_r/r·v_r/r              ← 简化后不需要
      + 2μ·∂u_z/∂r·∂v_r/∂z] dΩ             ← dw_laplace(πμr, v_r, u_z) [注意交叉项]

注意：u_r/r·v_r/r 项可以改写为：
∫ 2πr·λ·(u_r/r)·(v_r/r) dΩ = ∫ (2πλ/r)·u_r·v_r dΩ  ← dw_volume_dot(2πλ/r, v_r, u_r)

轴向方程（对 v_z 的系数）：

∫ 2πr [(λ+2μ)·∂u_z/∂z·∂v_z/∂z              ← dw_laplace(2πr(λ+2μ), v_z, u_z)
      + λ·∂u_r/∂r·∂v_z/∂z                   ← dw_laplace(2πrλ, v_z, u_r)
      + λ·u_r/r·∂v_z/∂z                     ← 合并到上面
      + 2μ·∂u_r/∂z·∂v_z/∂r] dΩ             ← dw_laplace(πμr, v_z, u_r)


7. SfePy Term 映射
-------------------
径向平衡：
    dw_laplace.i.Omega(m_rr.K, v_r, u_r)         K = 2πr(λ+2μ)
  + dw_laplace.i.Omega(m_prr.K, v_r, u_z)        K = 2πrλ
  + dw_volume_dot.i.Omega(m_hoop.c, v_r, u_r)   c = 2πλ/r
  + dw_laplace.i.Omega(m_shear.K, v_r, u_z)     K = πμr
  = 0

轴向平衡：
    dw_laplace.i.Omega(m_zz.K, v_z, u_z)         K = 2πr(λ+2μ)
  + dw_laplace.i.Omega(m_pzz.K, v_z, u_r)        K = 2πrλ
  + dw_laplace.i.Omega(m_shear.K, v_z, u_r)     K = πμr
  = 0

================================================================================

Key insight:
-----------
Rather than implementing new C code, we decompose the weak form into existing 
SfePy terms with spatially-varying coefficients that capture the 2πr weight 
and 1/r terms.

We use TWO SCALAR FIELDS (u_r and u_z) instead of a vector field, allowing us to:
- Apply different coefficients to each component
- Mix gradient terms (dw_laplace) with volume terms (dw_volume_dot)
- Handle the hoop strain u_r/r naturally

Test case:
----------
Domain: Annulus with r ∈ [0.1, 0.2], z ∈ [0, 0.3]
Manufactured solution: u_r = a*r, u_z = b*z
Expected: Machine precision for linear elements

Usage:
------
python3 -m sfepy.scripts.simple sfepy/examples/linear_elasticity/elastic_axisym_composition_v2.py
"""
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
