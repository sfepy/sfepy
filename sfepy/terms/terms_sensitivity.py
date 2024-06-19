import numpy as nm
from sfepy.linalg import dot_sequences
from sfepy.terms.terms_multilinear import ETermBase, sym2nonsym


def get_nonsym_grad_op(sgrad):
    nel, nqp, dim, _ = sgrad.shape
    grad_op = nm.zeros((nel, nqp, dim**2, dim**2), dtype=nm.float64)
    if dim == 3:
        grad_op[..., 0:3, 0:3] = sgrad
        grad_op[..., 3:6, 3:6] = sgrad
        grad_op[..., 6:9, 6:9] = sgrad
    elif dim == 2:
        grad_op[..., 0:2, 0:2] = sgrad
        grad_op[..., 2:4, 2:4] = sgrad
    else:
        grad_op = sgrad

    return grad_op


class ESDLinearElasticTerm(ETermBase):
    r"""
    Sensitivity analysis of the linear elastic term.
    :math:`D_{ijkl}` can be given in symmetric or non-symmetric form.

    :Definition:

    .. math::
        \int_{\Omega} \hat{D}_{ijkl}
        {\partial v_i \over \partial x_j}
        {\partial u_k \over \partial x_l}

    .. math::
        \hat{D}_{ijkl} = D_{ijkl}(\nabla \cdot \ul{\Vcal})
        - D_{ijkq}{\partial \Vcal_l \over \partial x_q}
        - D_{iqkl}{\partial \Vcal_j \over \partial x_q}

    :Arguments 1:
        - material : :math:`\ull{D}`
        - virtual/parameter_v : :math:`\ul{v}`
        - state/parameter_s : :math:`\ul{u}`
        - parameter_mv : :math:`\ul{\Vcal}`
    """
    name = 'de_sd_lin_elastic'
    arg_types = (('material', 'virtual', 'state', 'parameter_mv'),
                 ('material', 'parameter_1', 'parameter_2', 'parameter_mv'))
    arg_shapes = [{'material': 'S, S',
                  'virtual': ('D', 'state'), 'state': 'D',
                  'parameter_1': 'D', 'parameter_2': 'D',
                  'parameter_mv': 'D'}, {'material': 'D2, D2'}]
    modes = ('weak', 'eval')
    geometries = ['2_3', '2_4', '3_4', '3_8']

    def get_function(self, mat, vvar, svar, par_mv,
                     mode=None, term_mode=None, diff_var=None, **kwargs):

        grad_mv = self.get(par_mv, 'grad')
        div_mv = nm.trace(grad_mv, axis1=2, axis2=3)[..., None, None]
        grad_op = get_nonsym_grad_op(grad_mv)

        dim = grad_mv.shape[-1]
        mat_ns = mat if mat.shape[-1] == dim**2 else sym2nonsym(mat, [2, 3])
        aux = dot_sequences(mat_ns, grad_op, mode='AB')
        mat_sd = mat_ns * div_mv - aux - aux.transpose((0, 1, 3, 2))

        fun = self.make_function(
            'IK,v(i.j)->I,v(k.l)->K', (mat_sd, 'material'), vvar, svar,
            mode=mode, diff_var=diff_var,
        )

        return fun


class ESDPiezoCouplingTerm(ETermBase):
    r"""
    Sensitivity (shape derivative) of the piezoelectric coupling term.

    :Definition:

    .. math::
        \int_{\Omega} \hat{g}_{kij}\ e_{ij}(\ul{v}) \nabla_k p \mbox{ , }
        \int_{\Omega} \hat{g}_{kij}\ e_{ij}(\ul{u}) \nabla_k q

    .. math::
        \hat{g}_{kij} = g_{kij}(\nabla \cdot \ul{\Vcal})
        - g_{kil}{\partial \Vcal_j \over \partial x_l}
        - g_{lij}{\partial \Vcal_k \over \partial x_l}

    :Arguments 1:
        - material    : :math:`g_{kij}`
        - virtual/parameter_v : :math:`\ul{v}`
        - state/parameter_s : :math:`p`
        - parameter_mv : :math:`\ul{\Vcal}`

    :Arguments 2:
        - material : :math:`g_{kij}`
        - state    : :math:`\ul{u}`
        - virtual  : :math:`q`
        - parameter_mv : :math:`\ul{\Vcal}`
    """
    name = 'de_sd_piezo_coupling'
    arg_types = (('material', 'virtual', 'state', 'parameter_mv'),
                 ('material', 'state', 'virtual', 'parameter_mv'),
                 ('material', 'parameter_v', 'parameter_s', 'parameter_mv'))
    arg_shapes = {'material' : 'D, S',
                  'virtual/grad' : ('D', None), 'state/grad' : 1,
                  'virtual/div' : (1, None), 'state/div' : 'D',
                  'parameter_v' : 'D', 'parameter_s' : 1,
                  'parameter_mv': 'D'}
    modes = ('grad', 'div', 'eval')
    geometries = ['2_3', '2_4', '3_4', '3_8']

    def get_function(self, mat, vvar, svar, par_mv,
                     mode=None, term_mode=None, diff_var=None, **kwargs):
        grad_mv = self.get(par_mv, 'grad')
        div_mv = nm.trace(grad_mv, axis1=2, axis2=3)[..., None, None]
        grad_op = get_nonsym_grad_op(grad_mv)

        mat_ns = sym2nonsym(mat, [3])
        mat_sd = mat_ns * div_mv - dot_sequences(mat_ns, grad_op, mode='AB')\
               - dot_sequences(grad_mv, mat_ns, mode='ATB')

        expr = 'Ik,0.k,v(i.j)->I' if mode == 'div' else 'kI,v(i.j)->I,0.k'

        fun = self.make_function(
            expr, (mat_sd, 'material'), vvar, svar,
            mode=mode, diff_var=diff_var,
        )

        return fun


class ESDDiffusionTerm(ETermBase):
    r"""
    Diffusion sensitivity analysis term.

    :Definition:

    .. math::
        \int_{\Omega} \hat{K}_{ij} \nabla_i q\, \nabla_j p

    .. math::
        \hat{K}_{ij} = K_{ij}\left(
            \delta_{ik}\delta_{jl} \nabla \cdot \ul{\Vcal}
          - \delta_{ik}{\partial \Vcal_j \over \partial x_l}
          - \delta_{jl}{\partial \Vcal_i \over \partial x_k}\right)

    :Arguments:
        - material: :math:`K_{ij}`
        - virtual/parameter_1: :math:`q`
        - state/parameter_2: :math:`p`
        - parameter_mv: :math:`\ul{\Vcal}`
    """
    name = 'de_sd_diffusion'
    arg_types = (('material', 'virtual', 'state', 'parameter_mv'),
                 ('material', 'parameter_1', 'parameter_2', 'parameter_mv'))
    arg_shapes = {'material' : 'D, D', 'virtual': (1, 'state'), 'state': 1,
                  'parameter_1': 1, 'parameter_2': 1, 'parameter_mv': 'D'}
    modes = ('weak', 'eval')

    def get_function(self, mat, vvar, svar, par_mv,
                     mode=None, term_mode=None, diff_var=None, **kwargs):
        grad_mv = self.get(par_mv, 'grad')
        div_mv = nm.trace(grad_mv, axis1=2, axis2=3)[..., None, None]

        aux = dot_sequences(mat, grad_mv, mode='AB')
        mat_sd = mat * div_mv - aux - aux.transpose((0, 1, 3, 2))

        fun = self.make_function(
            'ij,0.i,0.j', (mat_sd, 'material'), vvar, svar,
            mode=mode, diff_var=diff_var,
        )

        return fun


class ESDStokesTerm(ETermBase):
    r"""
    Stokes problem coupling term. Corresponds to weak forms of gradient and
    divergence terms.

    :Definition:

    .. math::
        \int_{\Omega} p\, \hat{I}_{ij} {\partial v_i \over \partial x_j} \mbox{ , }
        \int_{\Omega} q\, \hat{I}_{ij} {\partial u_i \over \partial x_j}

    .. math::
        \hat{I}_{ij} = \delta_{ij} \nabla \cdot \Vcal
          - {\partial \Vcal_j \over \partial x_i}

    :Arguments 1:
        - virtual/parameter_v: :math:`\ul{v}`
        - state/parameter_s: :math:`p`
        - parameter_mv: :math:`\ul{\Vcal}`

    :Arguments 2:
        - state    : :math:`\ul{u}`
        - virtual  : :math:`q`
        - parameter_mv: :math:`\ul{\Vcal}`
    """
    name = 'de_sd_stokes'
    arg_types = (('opt_material', 'virtual', 'state', 'parameter_mv'),
                 ('opt_material', 'state', 'virtual', 'parameter_mv'),
                 ('opt_material', 'parameter_v', 'parameter_s', 'parameter_mv'))
    arg_shapes = [{'opt_material': '1, 1',
                   'virtual/grad': ('D', None), 'state/grad': 1,
                   'virtual/div': (1, None), 'state/div': 'D',
                   'parameter_v': 'D', 'parameter_s': 1, 'parameter_mv': 'D'},
                  {'opt_material': None}]
    modes = ('grad', 'div', 'eval')
    texpr = 'ij,i.j,0'

    def get_function(self, coef, vvar, svar, par_mv,
                     mode=None, term_mode=None, diff_var=None, **kwargs):
        grad_mv = self.get(par_mv, 'grad')
        div_mv = nm.trace(grad_mv, axis1=2, axis2=3)[..., None, None]

        mul = coef if coef is not None else 1.

        dim = grad_mv.shape[-2]
        mat_sd = (nm.eye(dim) * div_mv - grad_mv) * mul

        fun = self.make_function(
            self.texpr, (mat_sd, 'material'), vvar, svar,
            mode=mode, diff_var=diff_var
        )

        return fun


class ESDDivGradTerm(ETermBase):
    r"""
     Sensitivity (shape derivative) of diffusion term `de_div_grad`.

    :Definition:

    .. math::
        \int_{\Omega} \hat{I} \nabla \ul{v} : \nabla \ul{u} \mbox{ , }
        \int_{\Omega} \nu \hat{I}  \nabla \ul{v} : \nabla \ul{u}

    .. math::
        \hat{I}_{ijkl} =
            \delta_{ik}\delta_{jl} \nabla \cdot \ul{\Vcal}
          - \delta_{ik}\delta_{js} {\partial \Vcal_l \over \partial x_s}
          - \delta_{is}\delta_{jl} {\partial \Vcal_k \over \partial x_s}

    :Arguments:
        - material: :math:`\nu` (viscosity, optional)
        - virtual/parameter_1: :math:`\ul{v}`
        - state/parameter_2: :math:`\ul{u}`
        - parameter_mv: :math:`\ul{\Vcal}`
    """
    name = 'de_sd_div_grad'
    arg_types = (('opt_material', 'virtual', 'state', 'parameter_mv'),
                 ('opt_material', 'parameter_1', 'parameter_2',
                  'parameter_mv'))
    arg_shapes = [{'opt_material': '1, 1', 'virtual': ('D', 'state'),
                   'state': 'D', 'parameter_1': 'D', 'parameter_2': 'D',
                   'parameter_mv': 'D'},
                  {'opt_material': None}]
    modes = ('weak', 'eval')

    def get_function(self, mat, vvar, svar, par_mv,
                     mode=None, term_mode=None, diff_var=None, **kwargs):
        grad_mv = self.get(par_mv, 'grad')
        div_mv = nm.trace(grad_mv, axis1=2, axis2=3)[..., None, None]
        grad_op = get_nonsym_grad_op(grad_mv)

        mat_sd = nm.eye(grad_mv.shape[-2]**2) * div_mv\
               - grad_op - grad_op.transpose((0, 1, 3, 2))

        if mat is not None:
            mat_sd *= mat

        fun = self.make_function(
            'IK,v(i.j)->I,v(k.l)->K', (mat_sd, 'material'), vvar, svar,
            mode=mode, diff_var=diff_var,
        )

        return fun


class ESDDotTerm(ETermBase):
    r"""
    Sensitivity (shape derivative) of dot product of scalars or vectors.

    :Definition:

    .. math::
        \int_\Omega q p (\nabla \cdot \ul{\Vcal}) \mbox{ , }
        \int_\Omega (\ul{v} \cdot \ul{u}) (\nabla \cdot \ul{\Vcal})\\
        \int_\Omega c q p (\nabla \cdot \ul{\Vcal}) \mbox{ , }
        \int_\Omega c (\ul{v} \cdot \ul{u}) (\nabla \cdot \ul{\Vcal})\\
        \int_\Omega \ul{v} \cdot (\ull{M}\, \ul{u}) (\nabla \cdot \ul{\Vcal})

    :Arguments:
        - material: :math:`c` or :math:`\ull{M}` (optional)
        - virtual/parameter_1: :math:`q` or :math:`\ul{v}`
        - state/parameter_2: :math:`p` or :math:`\ul{u}`
        - parameter_mv : :math:`\ul{\Vcal}`
    """
    name = 'de_sd_dot'
    arg_types = (('opt_material', 'virtual', 'state', 'parameter_mv'),
                 ('opt_material', 'parameter_1', 'parameter_2',
                  'parameter_mv'))
    arg_shapes = [{'opt_material': '1, 1', 'virtual': (1, 'state'),
                   'state': 1, 'parameter_1': 1, 'parameter_2': 1,
                   'parameter_mv': 'D'},
                  {'opt_material': None},
                  {'opt_material': '1, 1', 'virtual': ('D', 'state'),
                   'state': 'D', 'parameter_1': 'D', 'parameter_2': 'D',
                   'parameter_mv': 'D'},
                  {'opt_material': 'D, D'},
                  {'opt_material': None}]
    modes = ('weak', 'eval')

    def get_function(self, mat, vvar, svar, par_mv, mode=None, term_mode=None,
                     diff_var=None, **kwargs):
        mat_sd = self.get(par_mv, 'div')

        if mat is not None:
            mat_sd = mat_sd * mat

        if mat_sd.shape[-1] > 1:
            fun = self.make_function(
                'ij,i,j', (mat_sd, 'material'), vvar, svar,
                mode=mode, diff_var=diff_var,
            )

        else:
            fun = self.make_function(
                '00,i,i', (mat_sd, 'material'), vvar, svar,
                mode=mode, diff_var=diff_var,
            )

        return fun


class ESDLinearTractionTerm(ETermBase):
    r"""
    Sensitivity of the linear traction term.

    :Definition:

    .. math::
        \int_{\Gamma} \ul{v} \cdot \left[\left(\ull{\hat{\sigma}}\,
        \nabla \cdot \ul{\cal{V}} - \ull{{\hat\sigma}}\, \nabla \ul{\cal{V}}
        \right)\ul{n}\right]

    .. math::
        \ull{\hat\sigma} = \ull{I} \mbox{ , }
        \ull{\hat\sigma} = c\,\ull{I} \mbox{ or }
        \ull{\hat\sigma} = \ull{\sigma}

    :Arguments:
        - material: :math:`c`, :math:`\ul{\sigma}`, :math:`\ull{\sigma}`
        - virtual/parameter: :math:`\ul{v}`
        - parameter_mv: :math:`\ul{\Vcal}`
    """
    name = 'de_sd_surface_ltr'
    arg_types = (('opt_material', 'virtual', 'parameter_mv'),
                 ('opt_material', 'parameter', 'parameter_mv'))
    arg_shapes = [{'opt_material': 'S, 1', 'virtual': ('D', None),
                   'parameter_mv': 'D', 'parameter': 'D'},
                  {'opt_material': None}, {'opt_material': '1, 1'},
                  {'opt_material': 'D, D'}]
    modes = ('weak', 'eval')
    integration = 'facet'

    def get_function(self, traction, vvar, par_mv, mode=None, term_mode=None,
                     diff_var=None, **kwargs):
        sg, _ = self.get_mapping(vvar)
        grad_mv = self.get(par_mv, 'grad', integration='facet_extra')
        div_mv = nm.trace(grad_mv, axis1=2, axis2=3)[..., None, None]

        _, n_qp, dim, _ = grad_mv.shape
        sym = (dim + 1) * dim // 2

        tdim, tdim2 = (None, None) if traction is None else traction.shape[2:]

        if tdim is None:
            trac = nm.tile(nm.eye(dim), (1, n_qp, 1, 1))
        elif tdim == 1:  # scalar value
            trac = nm.tile(nm.eye(dim), (1, n_qp, 1, 1)) * traction
        elif tdim == dim and tdim2 == dim:  # traction tensor - full
            trac = traction
        elif tdim == dim**2 and tdim2 == 1:  # traction tensor - full, vector
            trac = traction.reshape((-1, n_qp, dim, dim))
        elif tdim == sym and tdim2 == 1:  # traction tensor - symetric
            trac = sym2nonsym(traction, [2]).reshape((-1, n_qp, dim, dim))
        else:
            raise NotImplementedError

        aux = trac * div_mv - dot_sequences(trac, grad_mv, 'AB')
        mat_sd = dot_sequences(aux, sg.normal, 'AB')[..., 0]

        fun = self.make_function(
            'i,i', (mat_sd, 'opt_material'), vvar, mode=mode, diff_var=diff_var,
        )

        return fun


class ESDVectorDotGradScalarTerm(ESDStokesTerm):
    r"""
    Sensitivity of volume dot product of a vector and a gradient of scalar.

    :Definition:

    .. math::
        \int_{\Omega} \hat{I}_{ij} {\partial p \over \partial x_j}\, v_i
        \mbox{ , }
        \int_{\Omega} \hat{I}_{ij} {\partial q \over \partial x_j}\, u_i

    .. math::
        \hat{I}_{ij} = \delta_{ij} \nabla \cdot \Vcal
          - {\partial \Vcal_j \over \partial x_i}

    :Arguments 1:
        - virtual/parameter_v: :math:`\ul{v}`
        - state/parameter_s: :math:`p`
        - parameter_mv: :math:`\ul{\Vcal}`

    :Arguments 2:
        - state    : :math:`\ul{u}`
        - virtual  : :math:`q`
        - parameter_mv: :math:`\ul{\Vcal}`
    """
    name = 'de_sd_v_dot_grad_s'
    texpr = 'ij,i,0.j'
