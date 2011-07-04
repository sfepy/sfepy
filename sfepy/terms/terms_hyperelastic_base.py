import numpy as nm

from sfepy.base.base import assert_, Struct
from sfepy.terms.terms import Term, terms

class CouplingVectorScalarHE(Struct):

    def set_data_shape(self, apv, aps):
        """Set element data shape and checks dimensions of approximations."""
        if self.dof_conn_type == 'surface':
            self.data_shape_v = apv.get_s_data_shape(self.integral,
                                                     self.region.name)
            self.data_shape_s = aps.get_s_data_shape(self.integral,
                                                     self.region.name)
        else:
            self.data_shape_v = apv.get_v_data_shape(self.integral)
            self.data_shape_s = aps.get_v_data_shape(self.integral)

        assert_(aps.dim == (1,))
        assert_(apv.dim == (self.data_shape_v[2],))

    def get_shape_div(self, diff_var, chunk_size):
        n_el, n_qp, dim, n_eps = self.data_shape_s

        if diff_var is None:
            return (chunk_size, 1, n_eps, 1), 0
        elif diff_var == self.get_arg_name( 'state_p' ):
            return (chunk_size, 1, n_eps, n_eps), 2
        elif diff_var == self.get_arg_name( 'state' ):
            n_epv = self.data_shape_v[3]
            return (chunk_size, 1, n_eps, dim * n_epv), 1
        else:
            raise StopIteration

    def get_shape_grad(self, diff_var, chunk_size):
        n_el, n_qp, dim, n_epv = self.data_shape_v

        if diff_var is None:
            return (chunk_size, 1, dim * n_epv, 1), 0
        elif diff_var == self.get_arg_name( 'state' ):
            return (chunk_size, 1, dim * n_epv, dim * n_epv), 1
        elif diff_var == self.get_arg_name( 'state_p' ):
            n_eps = self.data_shape_s[3]
            return (chunk_size, 1, dim * n_epv, n_eps), 2
        else:
            raise StopIteration

class HyperElasticBase(Term):
    """
    Base class for all hyperelastic terms in TL/UL formulation.

    Note
    ----
    This is not a proper Term!

    `HyperElasticBase.__call__()` computes element contributions given either
    stress (-> rezidual) or tangent modulus (-> tangent sitffnes matrix),
    i.e. constitutive relation type (CRT) related data. The CRT data are
    computed in subclasses implementing particular CRT (e.g. neo-Hookean
    material), in self.compute_crt_data().
    Mode: 0 - total formulation, 1 - updated formulation
    """
    def __init__(self, *args, **kwargs):
        Term.__init__(self, *args, **kwargs)

        igs = self.region.igs
        self.stress_cache = {}.fromkeys(igs, None)

    def get_family_data(self, state, data_names):
        """
        Note
        ----
        `data_names` argument is ignored for now.
        """
        cache = state.evaluate_cache.setdefault('tl_common', {})

        vg, _, key = self.get_mapping(state, return_key=True)

        name = state.name
        data_key = key + (self.arg_steps[name], self.arg_derivatives[name])

        if data_key in cache:
            out = cache[data_key]

        else:
            out = self.compute_family_data(state)
            cache[data_key] = out

        return out

    def _call_hmode( self, diff_var = None, chunk_size = None, **kwargs ):
        term_mode, = self.get_kwargs( ['term_mode'], **kwargs )
        virtual, state, state_u = self.get_args( ['virtual', 'state', 'state_u'], **kwargs )
        ap, vg = self.get_approximation(virtual)

        self.set_data_shape(ap)
        shape, mode = self.get_shape(diff_var, chunk_size)

        cache = self.get_cache( self.function['finite_strain'][self.mode_ul], 0 )
        family_data = cache(self.family_data_names, self, 0, state=state_u)

        ig = self.char_fun.ig
        if term_mode is None:
            stress = self.crt_data.stress[ig]
            if stress is None:
                stress = self.compute_crt_data(family_data, 0, **kwargs)
                self.crt_data.stress[ig] = stress
            tan_mod = self.crt_data.tan_mod[ig]
            if tan_mod is None:
                tan_mod = self.compute_crt_data(family_data, 1, **kwargs)
                self.crt_data.tan_mod[ig] = tan_mod

            fun = self.function['element_contribution']
            mtxF, detF = cache(['F', 'detF'], self, 0, state=state_u)

            if mode == 0:
                vec = self.get_vector(state)
                for out, chunk in self.char_fun( chunk_size, shape ):
                    out2 = nm.zeros(out.shape[:-1] + (out.shape[-2],),
                                    dtype=nm.float64)
                    status1 = fun( out2, stress, tan_mod,
                                  mtxF, detF, vg, chunk, 1, self.mode_ul )
                    status2 = terms.he_residuum_from_mtx( out, out2, vec, ap.econn, chunk )
                    yield out, chunk, status1 or status2
            else:
                for out, chunk in self.char_fun( chunk_size, shape ):
                    status = fun( out, stress, tan_mod,
                                  mtxF, detF, vg, chunk, 1, self.mode_ul )
                    yield out, chunk, status

        elif term_mode == 'd_eval':
            raise NotImplementedError

    def _call_emode( self, diff_var = None, chunk_size = None, **kwargs ):
        term_mode, = self.get_kwargs( ['term_mode'], **kwargs )
        par1, par2, state_u = self.get_args(['parameter_1', 'parameter_2', 'state_u'], **kwargs)
        ap, vg = self.get_approximation(par1)

        self.set_data_shape(ap)
        n_el, n_qp, dim, n_ep = self.data_shape
        shape0 = (1, dim * n_ep, dim * n_ep)
        shape = (chunk_size, 1, 1, 1)

        cache = self.get_cache( self.function['finite_strain'][self.mode_ul], 0 )
        family_data = cache(self.family_data_names, self, 0, state=state_u)

        ig = self.char_fun.ig

        stress = self.crt_data.stress[ig]
        if stress is None:
            stress = self.compute_crt_data(family_data, 0, **kwargs)
            self.crt_data.stress[ig] = stress
        tan_mod = self.crt_data.tan_mod[ig]
        if tan_mod is None:
            tan_mod = self.compute_crt_data(family_data, 1, **kwargs)
            self.crt_data.tan_mod[ig] = tan_mod

        fun = self.function['element_contribution']
        mtxF, detF = cache(['F', 'detF'], self, 0, state=state_u)

        p1 = self.get_vector(par1)
        p2 = self.get_vector(par2)
        for out, chunk in self.char_fun( chunk_size, shape ):
            out2 = nm.zeros((out.shape[0],) + shape0, dtype=nm.float64)
            status1 = fun( out2, stress, tan_mod,
                           mtxF, detF, vg, chunk, 1, self.mode_ul )
            status2 = terms.he_eval_from_mtx(out, out2, p1, p2, ap.econn, chunk)
            out0 = nm.sum(out)

            yield out0, chunk, status1 or status2

    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        if self.call_mode is 0:
            return self._call_smode( diff_var, chunk_size, **kwargs )
        elif self.call_mode is 1:
            return self._call_hmode( diff_var, chunk_size, **kwargs )
        elif self.call_mode is 2:
            return self._call_emode( diff_var, chunk_size, **kwargs )
        else:
            raise NotImplementedError

class DeformationGradientTerm(Term):
    r"""
    :Description:
    Deformation gradient :math:`F` in quadrature points for
    `term_mode='def_grad'` (default) or the jacobian :math:`J` if
    `term_mode='jacobian'`.

    :Definition:
    .. math::
        \ull{F} = \pdiff{\ul{x}}{\ul{X}}|_{qp}
        = \ull{I} + \pdiff{\ul{u}}{\ul{X}}|_{qp} \;, \\
        \ul{x} = \ul{X} + \ul{u} \;, J = \det{(\ull{F})}

    :Arguments:
        state : :math:`\ul{u}`
    """
    name = 'dq_def_grad'
    arg_types = ('state',)

    function = staticmethod(terms.dq_def_grad)

    def __call__(self, diff_var=None, chunk_size=None, **kwargs):
        state, = self.get_args(**kwargs)
        term_mode = kwargs.get('term_mode', 'dq_def_grad')

        ap, vg = self.get_approximation(state)
        n_el, n_qp, dim, n_ep = ap.get_v_data_shape(self.integral)

        if diff_var is None:
            if term_mode == 'def_grad':
                shape = (chunk_size, n_qp, dim, dim)
                mode = 0

            elif term_mode == 'jacobian':
                shape = (chunk_size, n_qp, 1, 1)
                mode = 1

        else:
            raise StopIteration

        vec = state()
        for out, chunk in self.char_fun(chunk_size, shape):
            status = self.function(out, vec, vg, ap.econn, chunk, mode)
            yield out, chunk, status
