import numpy as nm

from sfepy.base.base import assert_
from sfepy.terms.terms import Term, terms, reorder_dofs_on_mirror
from sfepy.terms.terms_base import ScalarScalar

class MassVectorTerm(Term):
    r"""
    :Description:
    Vector field mass matrix/rezidual.

    :Definition:
    .. math::
        \int_{\Omega} \rho\ \ul{v} \cdot \ul{u}

    :Arguments:
        material : :math:`\rho`,
        virtual  : :math:`\ul{v}`,
        state    : :math:`\ul{u}`
    """
    name = 'dw_mass_vector'
    arg_types = ('material', 'virtual', 'state')

    function = staticmethod(terms.dw_mass)

    def get_fargs(self, mat, virtual, state,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        vg, _ = self.get_mapping(state)

        aux = nm.array([0], ndmin=4, dtype=nm.float64)
        if diff_var is None:
            val_qp = self.get(state, 'val')
            fmode = 0

        else:
            val_qp = aux
            fmode = 1

        return mat, val_qp, vg.bf, vg, fmode

class MassScalarTerm(Term):
    r"""
    :Description:
    Scalar field mass matrix/rezidual.

    :Definition:
    .. math::
        \int_{\Omega} q p

    :Arguments 1:
        virtual : :math:`q`,
        state   : :math:`p`

    :Arguments 2:
        parameter_1 : :math:`r`,
        parameter_2 : :math:`p`
    """
    name = 'dw_mass_scalar'
    arg_types = (('virtual', 'state'),
                 ('parameter_1', 'parameter_2'))
    modes = ('weak', 'eval')

    def check_shapes(self, virtual, state):
        assert_(virtual.n_components == 1)
        assert_(state.n_components == 1)

    def get_fargs(self, virtual, state,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        vg, _ = self.get_mapping(state)

        if mode == 'weak':
            aux = nm.array([0], ndmin=4, dtype=nm.float64)
            if diff_var is None:
                val_qp = self.get(state, 'val')
                fmode = 0

            else:
                val_qp = aux
                fmode = 1

            coef = kwargs.get('material')
            if coef is None:
                coef = nm.ones((1, val_qp.shape[1], 1, 1), dtype=nm.float64)

            return coef, val_qp, vg.bf, vg, fmode

        elif mode == 'eval':
            val_qp1 = self.get(virtual, 'val')
            val_qp2 = self.get(state, 'val')

            coef = kwargs.get('material')
            if coef is None:
                coef = nm.ones((1, val_qp.shape[1], 1, 1), dtype=nm.float64)

            return coef, val_qp1, val_qp2, vg.bf, vg

        else:
            raise ValueError('unsupported evaluation mode in %s! (%s)'
                             % (self.name, mode))

    def set_arg_types( self ):
        if self.mode == 'weak':
            self.function = terms.dw_mass_scalar

        else:
            self.function = terms.d_mass_scalar

class MassScalarWTerm(MassScalarTerm):
    r"""
    :Description:
    Scalar field mass matrix/rezidual weighted by a scalar function :math:`c`.

    :Definition:
    .. math::
        \int_{\Omega} c q p

    :Arguments 1:
        material : :math:`c`,
        virtual  : :math:`q`,
        state    : :math:`p`

    :Arguments 2:
        material    : :math:`c`,
        parameter_1 : :math:`r`,
        parameter_2 : :math:`p`
    """
    name = 'dw_mass_scalar_w'
    arg_types = (('material', 'virtual', 'state'),
                 ('material', 'parameter_1', 'parameter_2'))

    def check_shapes(self, mat, virtual, state):
        MassScalarTerm.check_shapes(self, virtual, state)

        n_el, n_qp, dim, n_en, n_c = self.get_data_shape(state)

        assert_(mat.shape[1:] == (n_qp, 1, 1))
#        assert_((mat.shape[0] == 1) or (mat.shape[0] == n_el))

    def get_fargs(self, material, virtual, state,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        fargs = MassScalarTerm.get_fargs(self, virtual, state,
                                         mode, term_mode, diff_var,
                                         material=material, **kwargs)
        return fargs

class MassScalarSurfaceTerm( ScalarScalar, Term ):
    r"""
    :Description:
    Scalar field mass matrix/rezidual on a surface.

    :Definition:
    .. math::
        \int_{\Gamma} q p

    :Arguments:
        virtual : :math:`q`,
        state   : :math:`p`
    """
    name = 'dw_surface_mass_scalar'
    arg_types = ('virtual', 'state')
    integration = 'surface'

    function = staticmethod(terms.dw_surf_mass_scalar)

    def get_fargs( self, diff_var = None, chunk_size = None, **kwargs ):
        virtual, state = self.get_args( ['virtual', 'state'], **kwargs )
        ap, sg = self.get_approximation(virtual)
        aps, sgs = self.get_approximation(state)

        self.set_data_shape( ap )
        shape, mode = self.get_shape( diff_var, chunk_size )

        vec = self.get_vector( state )
        if self.region.name in ap.surface_data:
            sd = ap.surface_data[self.region.name]
        else:
            sd = aps.surface_data[self.region.name]

        bf = ap.get_base( sd.face_type, 0, self.integral )

        is_trace = self.arg_traces[state.name]
        if is_trace:
            mirror_region, _, _ = self.region.get_mirror_region()
            sds = aps.surface_data[mirror_region.name]
            econn = sds.get_connectivity(state.is_surface)
            dc_type = self.get_dof_conn_type()
            ig = self.region.igs[0]
            rgnt = virtual.get_global_node_tab(dc_type, ig);
            cgnt = state.get_global_node_tab(dc_type, ig,
                                             is_trace=is_trace);
            econn = reorder_dofs_on_mirror(econn, cgnt, rgnt)
        else:
            econn = sd.get_connectivity(state.is_surface)

        if 'material' in self.arg_types:
            coef, = self.get_args(['material'], **kwargs)
        else:
            coef = nm.ones((1, self.data_shape[1], 1, 1), dtype=nm.float64)

        if state.is_real():
            fargs = coef, vec, 0, bf, sgs, econn
        else:
            ac = nm.ascontiguousarray
            fargs = [(coef, ac(vec.real), 0, bf, sgs, econn),
                     (coef, ac(vec.imag), 0, bf, sgs, econn)]
            mode += 1j

        return fargs, shape, mode

class MassScalarSurfaceWTerm(MassScalarSurfaceTerm):
    r"""
    :Description:
    Scalar field mass matrix/rezidual on a surface weighted by a scalar function.

    :Definition:
    .. math::
        \int_{\Gamma} c q p

    :Arguments:
        material : :math:`c`,
        virtual  : :math:`q`,
        state    : :math:`p`
    """
    name = 'dw_surface_mass_scalar_w'
    arg_types = ('material', 'virtual', 'state')

class BCNewtonTerm(MassScalarSurfaceTerm):
    r"""
    :Description:
    Newton boundary condition term.

    :Definition:
    .. math::
        \int_{\Gamma} \alpha q (p - p_{\rm outer})

    :Arguments:
        material_1 : :math:`\alpha`,
        material_2 : :math:`p_{\rm outer}`,
        virtual    : :math:`q`,
        state      : :math:`p`
    """
    name = 'dw_bc_newton'
    arg_types = ('material_1', 'material_2', 'virtual', 'state')

    def get_fargs( self, diff_var = None, chunk_size = None, **kwargs ):
        shift, = self.get_args(['material_2'], **kwargs)
        call = MassScalarSurfaceTerm.get_fargs
        fargs, shape, mode = call(self, diff_var, chunk_size, **kwargs)

        if nm.isreal(mode):
            fargs = (fargs[0] - shift,) + fargs[1:]
        else:
            raise NotImplementedError
        
        return fargs, shape, mode

    def __call__(self, diff_var=None, chunk_size=None, **kwargs):
        coef, = self.get_args(['material_1'], **kwargs)

        call = MassScalarSurfaceTerm.__call__
        for out, chunk, status in call(self, diff_var, chunk_size, **kwargs):
            out = coef * out
            yield out, chunk, status

class MassScalarFineCoarseTerm( Term ):
    r"""
    :Description:
    Scalar field mass matrix/rezidual for coarse to fine grid
    interpolation. Field :math:`p_H` belong to the coarse grid, test field
    :math:`q_h` to the fine grid.

    :Definition:
    .. math::
        \int_{\Omega} q_h p_H

    :Arguments:
        virtual : :math:`q_h`,
        state   : :math:`p_H`,
        iemaps  : coarse-fine element maps,
        pbase   : coarse base functions
    """
    name = 'dw_mass_scalar_fine_coarse'
    arg_types = ('virtual', 'state', 'iemaps', 'pbase' )

    function = staticmethod(terms.dw_mass_scalar_fine_coarse)
        
    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        virtual, state, iemaps, pbase = self.get_args( **kwargs )
        apr, vgr = virtual.get_current_approximation()
        apc, vgc = virtual.get_current_approximation()
        n_el, n_qp, dim, n_epr = apr.get_v_data_shape()
        
        if diff_var is None:
            shape = (chunk_size, 1, n_epr, 1)
            mode = 0
        elif diff_var == self.get_arg_name( 'state' ):
            n_epc = apc.get_v_data_shape()[3]
            shape = (chunk_size, 1, n_epr, n_epc)
            mode = 1
        else:
            raise StopIteration

        vec = state()

        cbfs = pbase[self.char_fun.ig]
        iemap = iemaps[self.char_fun.ig]
        for out, chunk in self.char_fun( chunk_size, shape ):
            status = self.function( out, vec, 0, apr.bf['v'], cbfs,
                                    vgr, apc.econn, iemap, chunk, mode )
            
            yield out, chunk, status
