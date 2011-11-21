import numpy as nm

from sfepy.base.base import use_method_with_name
from sfepy.terms.terms import Term, terms
from sfepy.terms.terms_base import ScalarScalar

class DiffusionSATerm(Term):
    r"""
    :Description:
    Diffusion sensitivity analysis term.

    :Definition:
    .. math::
        \int_{\Omega} \left[ (\dvg \ul{\Vcal}) K_{ij} \nabla_i q\, \nabla_j p - K_{ij} (\nabla_j \ul{\Vcal} \nabla q) \nabla_i p - K_{ij} \nabla_j q (\nabla_i \ul{\Vcal} \nabla p)\right]

    :Arguments:
        material : :math:`K_{ij}`,
        parameter_q: :math:`q`,
        parameter_p: :math:`p`,
        parameter_v: :math:`\ul{\Vcal}`,
    """
    name = 'd_diffusion_sa'
    arg_types = ('material', 'parameter_q', 'parameter_p', 'parameter_v')

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

class SurfaceLaplaceLayerTerm(ScalarScalar, Term):
    r"""
    :Description:
    Acoustic 'layer' term - derivatives in surface directions.

    :Definition:
    .. math::
        \int_{\Gamma} c \partial_\alpha \ul{q}\,\partial_\alpha \ul{p}, \alpha = 1,\dots,N-1

    :Arguments 1:
        material: :math:`c`,
        virtual:  :math:`q`,
        state:    :math:`p`

    :Arguments 2:
        material: :math:`c`,
        parameter_1 : :math:`q`,
        parameter_2 : :math:`p`
    """
    name = 'dw_surface_laplace'
    arg_types = [('material', 'virtual', 'state'),
                 ('material', 'parameter_2', 'parameter_1')]
    modes = ('weak', 'eval')
    integration = 'surface'

    functions = {'weak': terms.dw_surf_laplace,
                 'eval': terms.d_surf_laplace}

    def get_fargs_weak(self, diff_var = None, chunk_size = None, **kwargs):
        coef, virtual, state = self.get_args(**kwargs)
        ap, sg = self.get_approximation(virtual)
        aps, sgs = self.get_approximation(state)

        self.set_data_shape(ap)
        shape, mode = self.get_shape(diff_var, chunk_size)

        vec = self.get_vector(state)
        sd = aps.surface_data[self.region.name]

        bfg = ap.get_base(sd.face_type, 1, self.integral)
        econn = sd.get_connectivity(state.is_surface)

        if state.is_real():
            fargs = vec, coef, bfg, sgs, econn
        else:
            ac = nm.ascontiguousarray
            fargs = [(ac(vec.real), coef, bfg, sgs, econn),
                     (ac(vec.imag), coef, bfg, sgs, econn)]
            mode += 1j

        return fargs, shape, mode

    def get_fargs_eval(self, diff_var = None, chunk_size = None, **kwargs):
        coef, par2, par1 = self.get_args(**kwargs)
        ap, sg = self.get_approximation(par2)
        aps, sgs = self.get_approximation(par1)

        self.set_data_shape(ap)
        shape, mode = self.get_shape(diff_var, chunk_size)

        sd = aps.surface_data[self.region.name]
        bfg = ap.get_base(sd.face_type, 1, self.integral)
        econn = sd.get_connectivity(par1.is_surface)

        if par1.is_real() and par2.is_real():
            fargs = (par1(), par2(), coef, bfg, sgs, econn)
        else:
            ac = nm.ascontiguousarray
            p1 = self.get_vector(par1)
            p2 = self.get_vector(par2)
            fargs = [(ac(p1.real), ac(p2.real), coef, bfg, sgs, econn),
                     (ac(p1.imag), ac(p2.imag), coef, bfg, sgs, econn),
                     (ac(p1.real), ac(p2.imag), coef, bfg, sgs, econn),
                     (ac(p1.imag), ac(p2.real), coef, bfg, sgs, econn)]
            mode += 1j

        return fargs, (chunk_size, 1, 1, 1), mode

    def set_arg_types(self):
        if self.mode == 'weak':
            self.function = self.functions['weak']
            use_method_with_name( self, self.get_fargs_weak, 'get_fargs' )
        else:
            self.function = self.functions['eval']
            use_method_with_name( self, self.get_fargs_eval, 'get_fargs' )

class SurfaceCoupleLayerTerm(ScalarScalar, Term):
    r"""
    :Description:
    Acoustic 'layer' term - derivatives in surface directions.

    :Definition:
    .. math::
        \int_{\Gamma} c q\,\partial_\alpha p,
        \int_{\Gamma} c \partial_\alpha p\, q,
        \int_{\Gamma} c \partial_\alpha r\, s,\alpha = 1,\dots,N-1

    :Arguments 1:
        material: :math:`c`,
        virtual:  :math:`q`,
        state:    :math:`p`

    :Arguments 2:
        material: :math:`c`,
        virtual:  :math:`q`,
        state:    :math:`p`

    :Arguments 3:
        material: :math:`c`,
        parameter_1 : :math:`s`,
        parameter_2 : :math:`r`

    """
    name = 'dw_surface_lcouple'
    arg_types = (('material', 'virtual', 'state'),
                 ('material', 'state', 'virtual'),
                 ('material', 'parameter_1', 'parameter_2'))
    modes = ('bv_ns', 'nv_bs', 'eval')
    integration = 'surface'
    functions = {'weak': terms.dw_surf_lcouple,
                 'eval': terms.d_surf_lcouple}

    def get_fargs_weak(self, diff_var = None, chunk_size = None, **kwargs):
        if self.mode == 'nv_bs':
            coef, state, virtual = self.get_args(**kwargs)
        else:
            coef, virtual, state = self.get_args(**kwargs)

        ap, sg = self.get_approximation(virtual)
        aps, sgs = self.get_approximation(state)

        self.set_data_shape(ap)
        shape, mode = self.get_shape(diff_var, chunk_size)

        vec = self.get_vector(state)
        sd = aps.surface_data[self.region.name]

        bf = ap.get_base(sd.face_type, 0, self.integral)
        bfg =  ap.get_base(sd.face_type, 1, self.integral)
        econn = sd.get_connectivity(state.is_surface)

        aux = coef.shape
        if self.mode == 'nv_bs':
            bf, bfg = bfg, bf
            ncoef = coef.reshape((aux[0], aux[1], 1, nm.max(aux[2:])))

        else:
            ncoef = coef.reshape((aux[0], aux[1], nm.max(aux[2:]), 1))

        if state.is_real():
            fargs = vec, ncoef, bf, bfg, sgs, econn

        else:
            ac = nm.ascontiguousarray
            fargs = [(ac(vec.real), ncoef, bf, bfg, sgs, econn),
                     (ac(vec.imag), ncoef, bf, bfg, sgs, econn)]
            mode += 1j

        return fargs, shape, mode

    def get_fargs_eval(self, diff_var = None, chunk_size = None, **kwargs):
        coef, par2, par1 = self.get_args(**kwargs)
        ap, sg = self.get_approximation(par2)
        aps, sgs = self.get_approximation(par1)

        self.set_data_shape(ap)
        shape, mode = self.get_shape(diff_var, chunk_size)

        sd = aps.surface_data[self.region.name]
        bf = ap.get_base(sd.face_type, 0, self.integral)
        bfg = ap.get_base(sd.face_type, 1, self.integral)
        econn = sd.get_connectivity(par1.is_surface)

        aux = coef.shape
        ncoef = coef.reshape((aux[0], aux[1], 1, nm.max(aux[2:])))

        if par1.is_real() and par2.is_real():
            fargs = (par1(), par2(), ncoef, bf, bfg, sgs, econn)
        else:
            ac = nm.ascontiguousarray
            p1 = self.get_vector(par1)
            p2 = self.get_vector(par2)
            fargs = [(ac(p1.real), ac(p2.real), ncoef, bf, bfg, sgs, econn),
                     (ac(p1.imag), ac(p2.imag), ncoef, bf, bfg, sgs, econn),
                     (ac(p1.real), ac(p2.imag), ncoef, bf, bfg, sgs, econn),
                     (ac(p1.imag), ac(p2.real), ncoef, bf, bfg, sgs, econn)]
            mode += 1j

        return fargs, (chunk_size, 1, 1, 1), mode

    def set_arg_types(self):
        if self.mode == 'bv_ns' or self.mode == 'nv_bs':
            self.function = self.functions['weak']
            use_method_with_name( self, self.get_fargs_weak, 'get_fargs' )
        else:
            self.function = self.functions['eval']
            use_method_with_name( self, self.get_fargs_eval, 'get_fargs' )
