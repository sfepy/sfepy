import numpy as nm

from sfepy.terms.terms import Term, terms

class DiffusionSATerm(Term):
    r"""
    Diffusion sensitivity analysis term.

    :Definition:

    .. math::
        \int_{\Omega} \left[ (\dvg \ul{\Vcal}) K_{ij} \nabla_i q\, \nabla_j p -
        K_{ij} (\nabla_j \ul{\Vcal} \nabla q) \nabla_i p - K_{ij} \nabla_j q
        (\nabla_i \ul{\Vcal} \nabla p)\right]

    :Arguments:
        - material:    :math:`K_{ij}`
        - parameter_q: :math:`q`
        - parameter_p: :math:`p`
        - parameter_v: :math:`\ul{\Vcal}`
    """
    name = 'd_diffusion_sa'
    arg_types = ('material', 'parameter_q', 'parameter_p', 'parameter_v')
    arg_shapes = {'material' : 'D, D',
                  'parameter_q' : 1, 'parameter_p' : 1, 'parameter_v' : 'D'}

    function = staticmethod(terms.d_diffusion_sa)

    def get_fargs(self, mat, parameter_q, parameter_p, parameter_v,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        vg, _ = self.get_mapping(parameter_p)

        grad_q = self.get(parameter_q, 'grad')
        grad_p = self.get(parameter_p, 'grad')
        grad_v = self.get(parameter_v, 'grad')
        div_v = self.get(parameter_v, 'div')

        return grad_q, grad_p, grad_v, div_v, mat, vg

    def get_eval_shape(self, mat, parameter_q, parameter_p, parameter_v,
                       mode=None, term_mode=None, diff_var=None, **kwargs):
        n_el, n_qp, dim, n_en, n_c = self.get_data_shape(parameter_q)

        return (n_el, 1, 1, 1), parameter_q.dtype

class SurfaceLaplaceLayerTerm(Term):
    r"""
    Acoustic 'layer' term - derivatives in surface directions.

    :Definition:

    .. math::
        \int_{\Gamma} c \partial_\alpha \ul{q}\,\partial_\alpha \ul{p}, \alpha
        = 1,\dots,N-1

    :Arguments 1:
        - material: :math:`c`
        - virtual:  :math:`q`
        - state:    :math:`p`

    :Arguments 2:
        - material:    :math:`c`
        - parameter_1: :math:`q`
        - parameter_2: :math:`p`
    """
    name = 'dw_surface_laplace'
    arg_types = [('material', 'virtual', 'state'),
                 ('material', 'parameter_2', 'parameter_1')]
    arg_shapes = {'material' : '1, 1', 'virtual' : (1, 'state'),
                  'state' : 1, 'parameter_1' : 1, 'parameter_2' : 1}
    modes = ('weak', 'eval')
    integration = 'surface'

    def get_fargs(self, mat, virtual, state,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        ap, sg = self.get_approximation(virtual)
        aps, sgs = self.get_approximation(state)

        sd = aps.surface_data[self.region.name]
        bfg = ap.get_base(sd.face_type, 1, self.integral)

        elvol = nm.sum(sgs.det, axis=1)
        vol = nm.tile(elvol[:,nm.newaxis], (1, sgs.det.shape[1], 1, 1))

        if mode == 'weak':
            if diff_var is None:
                vec = self.get(state, 'val', bfun=bfg) / vol
                fmode = 0

            else:
                vec = nm.array([0], ndmin=4, dtype=nm.float64)
                fmode = 1

            return vec, mat, bfg, sgs, fmode

        elif mode == 'eval':
            vec1 = self.get(state, 'val', bfun=bfg) / vol
            vec2 = self.get(virtual, 'val', bfun=bfg) / vol

            return vec1, vec2, mat, sgs

    def get_eval_shape(self, mat, virtual, state,
                       mode=None, term_mode=None, diff_var=None, **kwargs):
        n_el, n_qp, dim, n_en, n_c = self.get_data_shape(state)

        return (n_el, 1, 1, 1), state.dtype

    def set_arg_types(self):

        if self.mode == 'weak':
            self.function = terms.dw_surf_laplace
        else:
            self.function = terms.d_surf_laplace

class SurfaceCoupleLayerTerm(Term):
    r"""
    Acoustic 'layer' term - derivatives in surface directions.

    :Definition:

    .. math::
        \int_{\Gamma} c q\,\partial_\alpha p,
        \int_{\Gamma} c \partial_\alpha p\, q,
        \int_{\Gamma} c \partial_\alpha r\, s,\alpha = 1,\dots,N-1

    :Arguments 1:
        - material: :math:`c`
        - virtual:  :math:`q`
        - state:    :math:`p`

    :Arguments 2:
        - material: :math:`c`
        - virtual:  :math:`q`
        - state:    :math:`p`

    :Arguments 3:
        - material:    :math:`c`
        - parameter_1: :math:`s`
        - parameter_2: :math:`r`
    """
    name = 'dw_surface_lcouple'
    arg_types = (('material', 'virtual', 'state'),
                 ('material', 'state', 'virtual'),
                 ('material', 'parameter_1', 'parameter_2'))
    arg_shapes = {'material' : '1, 1', 'virtual' : (1, 'state'),
                  'state' : 1, 'parameter_1' : 1, 'parameter_2' : 1}
    modes = ('bv_ns', 'nv_bs', 'eval')
    integration = 'surface'

    def get_fargs(self, mat, virtual, state,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        if self.mode == 'nv_bs':
            state, virtual = virtual, state

        ap, sg = self.get_approximation(virtual)
        aps, sgs = self.get_approximation(state)

        sd = aps.surface_data[self.region.name]
        bf = ap.get_base(sd.face_type, 0, self.integral)
        bfg = ap.get_base(sd.face_type, 1, self.integral)

        elvol = nm.sum(sgs.det, axis=1)
        vol = nm.tile(elvol[:,nm.newaxis], (1, sgs.det.shape[1], 1, 1))

        if mode == 'weak':
            if self.mode == 'nv_bs':
                bf, bfg = bfg, bf
                nmat = mat.reshape((mat.shape[0], mat.shape[1],
                                    1, nm.max(mat.shape[2:])))

            else:
                nmat = mat.reshape((mat.shape[0], mat.shape[1],
                                    nm.max(mat.shape[2:]), 1))

            if diff_var is None:
                if self.mode == 'nv_bs':
                    vec = self.get(state, 'val', bfun=bfg) / vol
                else:
                    vec = self.get(state, 'val')
                fmode = 0

            else:
                vec = nm.array([0], ndmin=4, dtype=nm.float64)
                fmode = 1

            return vec, nmat, bf, bfg, sgs, fmode

        elif mode == 'eval':
            vec1 = self.get(state, 'val')
            vec2 = self.get(virtual, 'val', bfun=bfg) / vol
            nmat = mat.reshape((mat.shape[0], mat.shape[1],
                                1, nm.max(mat.shape[2:])))

            return vec1, vec2, nmat, sgs

    def get_eval_shape(self, mat, virtual, state,
                       mode=None, term_mode=None, diff_var=None, **kwargs):
        n_el, n_qp, dim, n_en, n_c = self.get_data_shape(state)

        return (n_el, 1, 1, 1), state.dtype

    def set_arg_types(self):
        if self.mode == 'bv_ns' or self.mode == 'nv_bs':
            self.function = terms.dw_surf_lcouple
        else:
            self.function = terms.d_surf_lcouple
