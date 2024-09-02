"""
Proof-of-concept JAX-based terms supporting automatic differentiation.
"""
from functools import partial

import numpy as np

try:
    from jax import config
    config.update("jax_enable_x64", True)
    import jax
    import jax.numpy as jnp

except ImportError:
    print('install jax and jaxlib to use terms_jax.py!')
    raise

from sfepy.terms import Term
from sfepy.mechanics.matcoefs import lame_from_youngpoisson

def get_state_per_cell(term, state):
    aux = state(step=term.arg_steps[state.name],
                derivative=term.arg_derivatives[state.name])
    vec = aux.reshape((-1, state.n_components))
    econn = state.field.get_econn(term.act_integration, term.region)
    # Transpose is required to have sfepy order (DBD).
    cstate = vec[econn].transpose((0, 2, 1))

    return cstate

@jax.jit
def get_strain(gu):
    return 0.5 * (gu + gu.T)

@jax.jit
def get_stress(lam, mu, gu):
    strain = get_strain(gu)
    eye = jnp.eye(gu.shape[-1])
    stress = lam * jnp.trace(strain) * eye + 2 * mu * strain
    return stress

@jax.jit
def ceval_elasticity_l(lam, mu, vbfg, ubfg, det, cu):
    gu = cu @ ubfg.transpose((0, 2, 1)) # This does DBD.
    stress = jax.vmap(get_stress, in_axes=[None, None, 0])(lam, mu, gu)
    vbfgd = vbfg * det
    val = jnp.sum(stress @ vbfgd, axis=0)
    return val

@partial(jax.jit, static_argnames=['plane'])
def ceval_elasticity_yp(young, poisson, plane, vbfg, ubfg, det, cu):
    gu = cu @ ubfg.transpose((0, 2, 1))
    lam, mu = lame_from_youngpoisson(young, poisson, plane=plane)
    stress = jax.vmap(get_stress, in_axes=[None, None, 0])(lam, mu, gu)
    vbfgd = vbfg * det
    val = jnp.sum(stress @ vbfgd, axis=0)
    return val

eval_elasticity_l = jax.jit(jax.vmap(ceval_elasticity_l,
                                     in_axes=[None, None, 0, 0, 0, 0]))
eval_jac_elasticity_l = jax.jit(jax.vmap(jax.jacobian(ceval_elasticity_l, -1),
                                         in_axes=[None, None, 0, 0, 0, 0]))
eval_lam_elasticity_l = jax.vmap(jax.jacobian(ceval_elasticity_l, 0),
                                 in_axes=[None, None, 0, 0, 0, 0])
eval_mu_elasticity_l = jax.vmap(jax.jacobian(ceval_elasticity_l, 1),
                                in_axes=[None, None, 0, 0, 0, 0])

eval_elasticity_yp = jax.jit(
    jax.vmap(ceval_elasticity_yp,
             in_axes=[None, None, None, 0, 0, 0, 0]),
    static_argnames=['plane'],
)
eval_jac_elasticity_yp = jax.jit(
    jax.vmap(jax.jacobian(ceval_elasticity_yp, -1),
             in_axes=[None, None, None, 0, 0, 0, 0]),
    static_argnames=['plane'],
)
eval_young_elasticity_yp = jax.vmap(jax.jacobian(ceval_elasticity_yp, 0),
                                    in_axes=[None, None, None, 0, 0, 0, 0])
eval_poisson_elasticity_yp = jax.vmap(jax.jacobian(ceval_elasticity_yp, 1),
                                      in_axes=[None, None, None, 0, 0, 0, 0])

class LinearElasticLADTerm(Term):
    r"""
    Homogeneous isotropic linear elasticity term differentiable w.r.t. material
    parameters :math:`\lambda`, :math:`\mu` (Lamé's parameters).

    :Definition:

    .. math::
        \int_{\Omega} D_{ijkl}\ e_{ij}(\ul{v}) e_{kl}(\ul{u})\\ \mbox{ with } \\
        D_{ijkl} = \mu (\delta_{ik} \delta_{jl}+\delta_{il} \delta_{jk}) +
        \lambda \ \delta_{ij} \delta_{kl}

    :Arguments:
        - material_1: :math:`\lambda` (Lamé's first parameter)
        - material_2: :math:`\mu` (Lamé's second parameter, shear modulus)
        - virtual/parameter_1: :math:`\ul{v}`
        - state/parameter_2: :math:`\ul{u}`
    """
    name = 'dw_lin_elastic_l_ad'
    arg_types = (('material_1', 'material_2', 'virtual', 'state'),)
    arg_shapes = {'material_1' : '1, 1', 'material_2' : '1, 1',
                  'virtual' : ('D', 'state'), 'state' : 'D'}
    modes = ('weak',)
    diff_info = {'material_1' : 1, 'material_2' : 1}

    def get_fargs(self, material1, material2, virtual, state,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        vgmap, _ = self.get_mapping(virtual)
        sgmap, _ = self.get_mapping(state)

        cu = get_state_per_cell(self, state)

        lam = material1[0, 0, 0, 0]
        mu = material2[0, 0, 0, 0]
        if diff_var is None:
            fun = eval_elasticity_l

        elif diff_var == state.name:
            fun = eval_jac_elasticity_l

        elif diff_var == 'material_1':
            fun = eval_lam_elasticity_l

        elif diff_var == 'material_2':
            fun = eval_mu_elasticity_l

        else:
            raise ValueError

        fargs = [lam, mu, vgmap.bfg, sgmap.bfg, vgmap.det, cu]

        return fun, fargs

    @staticmethod
    def function(out, fun, fargs):
        out[:] = np.asarray(fun(*fargs).reshape(out.shape))
        return 0

class LinearElasticYPADTerm(Term):
    r"""
    Homogeneous isotropic linear elasticity term differentiable w.r.t. material
    parameters :math:`E` (Young's modulus), :math:`\nu` (Poisson's ratio).

    :Definition:

    .. math::
        \int_{\Omega} D_{ijkl}\ e_{ij}(\ul{v}) e_{kl}(\ul{u})\\ \mbox{ with } \\
        D_{ijkl} = \mu (\delta_{ik} \delta_{jl}+\delta_{il} \delta_{jk}) +
        \lambda \ \delta_{ij} \delta_{kl}, \\ \mbox{ where } \\
        \lambda = E \nu / ((1 + \nu)(1 - 2\nu)), \\ \mu = E / 2(1 + \nu)

    :Arguments:
        - material_1: :math:`E`
        - material_2: :math:`\nu`
        - virtual/parameter_1: :math:`\ul{v}`
        - state/parameter_2: :math:`\ul{u}`
    """
    name = 'dw_lin_elastic_yp_ad'
    arg_types = (('material_1', 'material_2', 'virtual', 'state'),)
    arg_shapes = {'material_1' : '1, 1', 'material_2' : '1, 1',
                  'virtual' : ('D', 'state'), 'state' : 'D'}
    modes = ('weak',)
    diff_info = {'material_1' : 1, 'material_2' : 1}

    def get_fargs(self, material1, material2, virtual, state,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        vgmap, _ = self.get_mapping(virtual)
        sgmap, _ = self.get_mapping(state)

        cu = get_state_per_cell(self, state)

        young = material1[0, 0, 0, 0]
        poisson = material2[0, 0, 0, 0]
        if diff_var is None:
            fun = eval_elasticity_yp

        elif diff_var == state.name:
            fun = eval_jac_elasticity_yp

        elif diff_var == 'material_1':
            fun = eval_young_elasticity_yp

        elif diff_var == 'material_2':
            fun = eval_poisson_elasticity_yp

        else:
            raise ValueError

        fargs = [young, poisson, 'strain', vgmap.bfg, sgmap.bfg, vgmap.det, cu]

        return fun, fargs

    @staticmethod
    def function(out, fun, fargs):
        out[:] = np.asarray(fun(*fargs).reshape(out.shape))
        return 0

@jax.jit
def ceval_mass(density, vbf, ubf, det, cu):
    ru = density * cu @ ubf.transpose((0, 2, 1))
    vbfd = vbf * det
    val = jnp.sum(ru @ vbfd, axis=0)
    return val

eval_mass = jax.jit(jax.vmap(ceval_mass,
                             in_axes=[None, 0, 0, 0, 0]))
eval_jac_mass = jax.jit(jax.vmap(jax.jacobian(ceval_mass, -1),
                                 in_axes=[None, 0, 0, 0, 0]))
eval_density_mass = jax.vmap(jax.jacobian(ceval_mass, 0),
                             in_axes=[None, 0, 0, 0, 0])

class MassADTerm(Term):
    r"""
    Homogeneous mass term differentiable w.r.t. the material parameter.

    :Definition:

    .. math::
        \int_{\cal{D}} \rho \ul{v} \cdot \ul{u}

    :Arguments:
        - material_1: :math:`\rho`
        - virtual: :math:`\ul{v}`
        - state: :math:`\ul{u}`
    """
    name = 'dw_mass_ad'
    arg_types = (('material', 'virtual', 'state'),)
    arg_shapes = [{'material' : '1, 1',
                   'virtual' : ('D', 'state'), 'state' : 'D'},
                  {'virtual' : (1, 'state'), 'state' : 1}]
    modes = ('weak',)
    diff_info = {'material' : 1}
    integration = ('cell', 'facet')

    def get_fargs(self, material_density, virtual, state,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        vgmap, _ = self.get_mapping(virtual)
        sgmap, _ = self.get_mapping(state)

        cu = get_state_per_cell(self, state)

        density = material_density[0, 0, 0, 0]
        if diff_var is None:
            fun = eval_mass

        elif diff_var == state.name:
            fun = eval_jac_mass

        elif diff_var == 'material':
            fun = eval_density_mass

        else:
            raise ValueError

        bs = partial(np.broadcast_to, shape=(cu.shape[0],) + vgmap.bf.shape[1:])
        fargs = [density, bs(vgmap.bf), bs(sgmap.bf), vgmap.det, cu]

        return fun, fargs

    @staticmethod
    def function(out, fun, fargs):
        out[:] = np.asarray(fun(*fargs).reshape(out.shape))
        return 0

def _get_Eq(gu):
    return 0.5 * (gu + gu.T + gu.T @ gu)

def _get_E(bfg, cu):
    gu = cu @ bfg.transpose((0, 2, 1)) # This does DBD.
    val = jax.vmap(_get_Eq, in_axes=[0])(gu)
    return val

_get_dEv = jax.jacobian(_get_E, -1)

def get_neohook_strain_energy(mu, C):
    # (*) 0.5 *
    dim = C.shape[-1]
    trc  = jnp.trace(C)
    if dim == 2:
        trc += 1.0 # Plane strain.

    W = mu * (jnp.linalg.det(C)**(-1.0/3.0) * trc - 3.0)
    return W

_get_neohook_stress_2pk = jax.grad(get_neohook_strain_energy, -1)

def get_neohook_stress_2pk(mu, gu):
    eye = np.eye(gu.shape[-1])
    F = gu + eye
    # (*) 2 *
    S = _get_neohook_stress_2pk(mu, F.T @ F)
    return S

# This is very slow, worsens with problem size w.r.t. dw_tl_neohook and order 2
# is out of memory.
# @jax.jit
def ceval_neohook0(mu, vbfg, ubfg, det, cu):
    gu = cu @ ubfg.transpose((0, 2, 1)) # This does DBD.
    stress = jax.vmap(get_neohook_stress_2pk, in_axes=[None, 0])(mu, gu)
    dEv = _get_dEv(vbfg, cu) * det[..., None, None]
    val = jnp.sum(stress[..., None, None] * dEv, axis=[0, 1, 2])
    return val

def get_neohook_strain_energy_f(mu, F):
    dim = F.shape[-1]
    trf = jnp.trace(F.T @ F)
    if dim == 2:
        trf += 1.0 # Plane strain

    W = 0.5 * mu * (jnp.linalg.det(F)**(-2.0/3.0) * trf - 3.0)
    return W

_get_neohook_stress_1pk = jax.grad(get_neohook_strain_energy_f, -1)

def get_neohook_stress_1pk(mu, gu):
    eye = np.eye(gu.shape[-1])
    F = gu + eye
    S = _get_neohook_stress_1pk(mu, F)
    return S

# This is slow, improves with problem size w.r.t. sfepy and order 2 works, with
# refine about 5x slower than sfepy.
@jax.jit
def ceval_neohook(mu, vbfg, ubfg, det, cu):
    gu = cu @ ubfg.transpose((0, 2, 1)) # This does DBD.
    stress = jax.vmap(get_neohook_stress_1pk, in_axes=[None, 0])(mu, gu)
    vbfgd = vbfg * det
    val = jnp.sum(stress @ vbfgd, axis=0)
    return val

eval_neohook = jax.jit(jax.vmap(ceval_neohook,
                                in_axes=[None, 0, 0, 0, 0]))
eval_jac_neohook = jax.jit(jax.vmap(jax.jacobian(ceval_neohook, -1),
                                    in_axes=[None, 0, 0, 0, 0]))
eval_mu_neohook = jax.jit(jax.vmap(jax.jacobian(ceval_neohook, 0),
                                   in_axes=[None, 0, 0, 0, 0]))

class NeoHookeanTLADTerm(Term):
    r"""
    Homogeneous Hyperelastic neo-Hookean term differentiable w.r.t. the
    material parameter. Effective stress :math:`S_{ij} = \mu
    J^{-\frac{2}{3}}(\delta_{ij} - \frac{1}{3}C_{kk}C_{ij}^{-1})`.

    :Definition:

    .. math::
        \int_{\Omega} S_{ij}(\ul{u}) \delta E_{ij}(\ul{u};\ul{v})

    :Arguments:
        - material : :math:`\mu`
        - virtual  : :math:`\ul{v}`
        - state    : :math:`\ul{u}`
    """
    name = 'dw_tl_he_neohook_ad'
    arg_types = ('material', 'virtual', 'state')
    arg_shapes = {'material' : '1, 1', 'virtual' : ('D', 'state'),
                  'state' : 'D'}
    modes = ('weak',)
    diff_info = {'material' : 1}
    geometries = ['2_3', '2_4', '3_4', '3_8']

    def get_fargs(self, material, virtual, state,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        vgmap, _ = self.get_mapping(virtual)
        sgmap, _ = self.get_mapping(state)

        cu = get_state_per_cell(self, state)

        mu = material[0, 0, 0, 0]
        if diff_var is None:
            fun = eval_neohook

        elif diff_var == state.name:
            fun = eval_jac_neohook

        elif diff_var == 'material':
            fun = eval_mu_neohook

        else:
            raise ValueError

        fargs = [mu, vgmap.bfg, sgmap.bfg, vgmap.det, cu]
        fargs = [jax.device_put(val) for val in fargs]
        return fun, fargs

    @staticmethod
    def function(out, fun, fargs):
        out[:] = np.asarray(fun(*fargs).reshape(out.shape))
        return 0

def get_ogden_strain_energy_f(mu, alpha, F):
    J = jnp.linalg.det(F)
    C = J**(-2.0/3.0) * F.T @ F # Isochoric C.
    lams2, Ns = jnp.linalg.eigh(C)
    a5 = 0.5 * alpha
    W = (mu / alpha) * (jnp.sum(lams2**a5) - 3.0)
    return W

_get_ogden_stress_1pk = jax.grad(get_ogden_strain_energy_f, -1)

def get_ogden_stress_1pk(mu, alpha, gu):
    eye = np.eye(gu.shape[-1])
    F = gu + eye
    S = _get_ogden_stress_1pk(mu, alpha, F)
    return S

@jax.jit
def ceval_ogden(mu, alpha, vbfg, ubfg, det, cu):
    gu = cu @ ubfg.transpose((0, 2, 1)) # This does DBD.
    stress = jax.vmap(get_ogden_stress_1pk,
                      in_axes=[None, None, 0])(mu, alpha, gu)
    vbfgd = vbfg * det
    val = jnp.sum(stress @ vbfgd, axis=0)
    return val

eval_ogden = (jax.vmap(ceval_ogden,
                              in_axes=[None, None, 0, 0, 0, 0]))
eval_jac_ogden = (jax.vmap(jax.jacobian(ceval_ogden, -1),
                                  in_axes=[None, None, 0, 0, 0, 0]))
eval_mu_ogden = jax.jit(jax.vmap(jax.jacobian(ceval_ogden, 0),
                                 in_axes=[None, None, 0, 0, 0, 0]))
eval_alpha_ogden = jax.jit(jax.vmap(jax.jacobian(ceval_ogden, 1),
                                    in_axes=[None, None, 0, 0, 0, 0]))

class OgdenTLADTerm(Term):
    r"""
    Homogeneous hyperelastic Ogden model term differentiable w.r.t. the
    material parameters, with the strain energy density

    .. math::
        W = \frac{\mu}{\alpha} \, \left(
            \bar\lambda_1^{\alpha} + \bar\lambda_2^{\alpha}
             + \bar\lambda_3^{\alpha} - 3 \right) \; ,

    where :math:`\lambda_k, k=1, 2, 3` are the principal stretches, whose
    squares are the principal values of the right Cauchy-Green deformation
    tensor :math:`\mathbf{C}`. For more details see :class:`OgdenTLTerm
    <sfepy.terms.terms_hyperelastic_tl.OgdenTLTerm>`.

    WARNING: The current implementation fails to compute the tangent matrix
    when :math:`\mathbf{C}` has multiple eigenvalues (e.g. zero deformation).
    In that case nans are returned, as a result of dividing by zero. See [1],
    Section 11.2.3, page 385.

    [1] Borst, R. de, Crisfield, M.A., Remmers, J.J.C., Verhoosel, C.V., 2012.
    Nonlinear Finite Element Analysis of Solids and Structures, 2nd edition.
    ed. Wiley, Hoboken, NJ.

    :Definition:

    .. math::
        \int_{\Omega} S_{ij}(\ul{u}) \delta E_{ij}(\ul{u};\ul{v})

    :Arguments:
        - material_1 : :math:`\mu`
        - material_2 : :math:`\alpha`
        - virtual    : :math:`\ul{v}`
        - state      : :math:`\ul{u}`
    """
    name = 'dw_tl_he_ogden_ad'
    arg_types = ('material_mu', 'material_alpha', 'virtual', 'state')
    arg_shapes = {'material_mu' : '1, 1', 'material_alpha' : '1, 1',
                  'virtual' : ('D', 'state'), 'state' : 'D'}
    modes = ('weak',)
    diff_info = {'material_mu' : 1, 'material_alpha' : 1}
    geometries = ['3_4', '3_8']

    def get_fargs(self, material_mu, material_alpha, virtual, state,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        vgmap, _ = self.get_mapping(virtual)
        sgmap, _ = self.get_mapping(state)

        cu = get_state_per_cell(self, state)

        mu = material_mu[0, 0, 0, 0]
        alpha = material_alpha[0, 0, 0, 0]

        if diff_var is None:
            fun = eval_ogden

        elif diff_var == state.name:
            fun = eval_jac_ogden

        elif diff_var == 'material_mu':
            fun = eval_mu_ogden

        elif diff_var == 'material_alpha':
            fun = eval_alpha_ogden

        else:
            raise ValueError

        fargs = [mu, alpha, vgmap.bfg, sgmap.bfg, vgmap.det, cu]
        fargs = [jax.device_put(val) for val in fargs]
        return fun, fargs

    @staticmethod
    def function(out, fun, fargs):
        out[:] = np.asarray(fun(*fargs).reshape(out.shape))
        return 0
