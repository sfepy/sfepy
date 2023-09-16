"""
Proof-of-concept JAX-based terms supporting automatic differentiation.
"""
from functools import partial

import numpy as np

try:
    from jax.config import config
    config.update("jax_enable_x64", True)
    import jax
    import jax.numpy as jnp

except ImportError:
    print('install jax and jaxlib to use terms_jax.py!')
    raise

from sfepy.terms import Term
from sfepy.mechanics.matcoefs import lame_from_youngpoisson

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
        vgmap, _ = self.get_mapping(state)
        sgmap, _ = self.get_mapping(state)

        vecu = jnp.array(state().reshape((-1, vgmap.dim)))
        econn = state.field.get_econn(self.integration, self.region)
        # Transpose is required to have sfepy order (DBD).
        cu = vecu[econn].transpose((0, 2, 1))

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
        vgmap, _ = self.get_mapping(state)
        sgmap, _ = self.get_mapping(state)

        vecu = jnp.array(state().reshape((-1, vgmap.dim)))
        econn = state.field.get_econn(self.integration, self.region)
        # Transpose is required to have sfepy order (DBD).
        cu = vecu[econn].transpose((0, 2, 1))

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
    arg_shapes = {'material' : '1, 1',
                  'virtual' : ('D', 'state'), 'state' : 'D'}
    modes = ('weak',)
    diff_info = {'material' : 1}

    def get_fargs(self, material_density, virtual, state,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        vgmap, _ = self.get_mapping(state)
        sgmap, _ = self.get_mapping(state)

        vecu = jnp.array(state().reshape((-1, vgmap.dim)))
        econn = state.field.get_econn(self.integration, self.region)
        # Transpose is required to have sfepy order (DBD).
        cu = vecu[econn].transpose((0, 2, 1))

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
