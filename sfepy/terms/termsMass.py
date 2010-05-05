from sfepy.terms.terms import *
from sfepy.terms.terms_base import VectorVector, ScalarScalar

class MassTerm( VectorVector, Term ):
    r"""
    :Description:
    Inertial forces term.

    :Definition:
    .. math::
        \int_{\Omega} \rho \ul{v} \cdot \frac{\ul{u} - \ul{u}_0}{\dt}

    :Arguments:
    material : :math:`\rho`,
    ts.dt : :math:`\dt`,
    parameter : :math:`\ul{u}_0`"""
    name = 'dw_mass'
    arg_types = ('ts', 'material', 'virtual', 'state', 'parameter')
    geometry = [(Volume, 'virtual')]

    def __init__(self, name, sign, **kwargs):
        Term.__init__(self, name, sign, function=terms.dw_mass, **kwargs)

    def get_fargs(self, diff_var=None, chunk_size=None, **kwargs):
        ts, mat, virtual, state, state0 = self.get_args(**kwargs)        
        ap, vg = virtual.get_approximation(self.get_current_group(), 'Volume')

        self.set_data_shape(ap)
        shape, mode = self.get_shape(diff_var, chunk_size)

        dvec = state() - state0()
        rhodt = mat / ts.dt
        bf = ap.get_base('v', 0, self.integral_name)

        fargs = (rhodt, dvec, 0, bf, vg, ap.econn)
        return fargs, shape, mode

class MassVectorTerm( MassTerm ):
    r"""
    :Description:
    Vector field mass matrix/rezidual.

    :Definition:
    .. math::
        \int_{\Omega} \rho\ \ul{v} \cdot \ul{u}
    """
    name = 'dw_mass_vector'
    arg_types = ('material', 'virtual', 'state')
    geometry = [(Volume, 'virtual')]

    def get_fargs(self, diff_var=None, chunk_size=None, **kwargs):
        mat, virtual, state = self.get_args(**kwargs)        
        ap, vg = virtual.get_approximation(self.get_current_group(), 'Volume')

        self.set_data_shape(ap)
        shape, mode = self.get_shape(diff_var, chunk_size)

        vec = self.get_vector(state)
        bf = ap.get_base('v', 0, self.integral_name)
        fargs = (mat, vec, 0, bf, vg, ap.econn)
        return fargs, shape, mode

class MassScalarTerm(ScalarScalar, Term):
    r"""
    :Description:
    Scalar field mass matrix/rezidual.

    :Definition:
    .. math::
        \int_{\Omega} q p
    """
    name = 'dw_mass_scalar'
    arg_types = (('virtual', 'state'),
                 ('parameter_1', 'parameter_2')) 
    geometry = ([(Volume, 'virtual')],
                [(Volume, 'parameter_1')])
    modes = ('weak', 'eval')
    functions = {'weak': terms.dw_mass_scalar,
                 'eval': terms.d_mass_scalar}

    def get_fargs_weak( self, diff_var = None, chunk_size = None, **kwargs ):
        virtual, state = self.get_args( ['virtual', 'state'], **kwargs )
        ap, vg = virtual.get_approximation( self.get_current_group(), 'Volume' )

        self.set_data_shape( ap )
        shape, mode = self.get_shape( diff_var, chunk_size )

        vec = self.get_vector( state )
        bf = ap.get_base( 'v', 0, self.integral_name )

        if 'material' in self.arg_types:
            coef, = self.get_args(['material'], **kwargs)

        else:
            coef = nm.ones((1, self.data_shape[1], 1, 1), dtype=nm.float64)

        if state.is_real():
            fargs = coef, vec, bf, vg, ap.econn
        else:
            ac = nm.ascontiguousarray
            fargs = [(coef, ac( vec.real ), bf, vg, ap.econn),
                     (coef, ac( vec.imag ), bf, vg, ap.econn)]
            mode += 1j
            
        return fargs, shape, mode

    def get_fargs_eval( self, diff_var = None, chunk_size = None, **kwargs ):
        par1, par2 = self.get_args(['parameter_1', 'parameter_2'], **kwargs)
        ap, vg = par1.get_approximation( self.get_current_group(), 'Volume' )
        self.set_data_shape( ap )
        bf = ap.get_base( 'v', 0, self.integral_name )

        if 'material' in self.arg_types:
            coef, = self.get_args(['material'], **kwargs)

        else:
            coef = nm.ones((1, self.data_shape[1], 1, 1), dtype=nm.float64)

        return (coef, par1(), par2(), bf, vg, ap.econn), (chunk_size, 1, 1, 1), 0

    def set_arg_types( self ):
        if self.mode == 'weak':
            self.function = self.functions['weak']
            use_method_with_name( self, self.get_fargs_weak, 'get_fargs' )
        else:
            self.function = self.functions['eval']
            use_method_with_name( self, self.get_fargs_eval, 'get_fargs' )

class MassScalarWTerm(MassScalarTerm):
    r"""
    :Description:
    Scalar field mass matrix/rezidual weighted by a scalar function :math:`c`.

    :Definition:
    .. math::
        \int_{\Omega} c q p
    """
    name = 'dw_mass_scalar_w'
    arg_types = (('material', 'virtual', 'state'),
                 ('material', 'parameter_1', 'parameter_2')) 
    geometry = ([(Volume, 'virtual')],
                [(Volume, 'parameter_1')])

class MassScalarSurfaceTerm( ScalarScalar, Term ):
    r"""
    :Description:
    Scalar field mass matrix/rezidual on a surface.

    :Definition:
    .. math::
        \int_{\Gamma} q p
    """
    name = 'dw_surface_mass_scalar'
    arg_types = ('virtual', 'state')
    geometry = [(Surface, 'virtual')]

    def __init__(self, name, sign, **kwargs):
        Term.__init__(self, name, sign, dof_conn_type='surface',
                      function=terms.dw_surf_mass_scalar, **kwargs)

    def get_fargs( self, diff_var = None, chunk_size = None, **kwargs ):
        virtual, state = self.get_args( ['virtual', 'state'], **kwargs )
        ap, sg = virtual.get_approximation( self.get_current_group(), 'Surface' )
        
        self.set_data_shape( ap )
        shape, mode = self.get_shape( diff_var, chunk_size )

        vec = self.get_vector( state )
        sd = ap.surface_data[self.region.name]
        bf = ap.get_base( sd.face_type, 0, self.integral_name )

        if 'material' in self.arg_types:
            coef, = self.get_args(['material'], **kwargs)

        else:
            coef = nm.ones((1, self.data_shape[1], 1, 1), dtype=nm.float64)

        if state.is_real():
            fargs = coef, vec, 0, bf, sg, sd.econn

        else:
            ac = nm.ascontiguousarray
            fargs = [(coef, ac(vec.real), 0, bf, sg, sd.econn),
                     (coef, ac(vec.imag), 0, bf, sg, sd.econn)]
            mode += 1j

        return fargs, shape, mode

class MassScalarSurfaceWTerm(MassScalarSurfaceTerm):
    r"""
    :Description:
    Scalar field mass matrix/rezidual on a surface weighted by a scalar function.

    :Definition:
    .. math::
        \int_{\Gamma} c q p
    """
    name = 'dw_surface_mass_scalar_w'
    arg_types = ('material', 'virtual', 'state')
    geometry = [(Surface, 'virtual')]

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
    virtual : :math:`q`,
    state : :math:`p`
    """
    name = 'dw_bc_newton'
    arg_types = ('material_1', 'material_2', 'virtual', 'state')
    geometry = [(Surface, 'virtual'), (Surface, 'state')]

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
    """
    name = 'dw_mass_scalar_fine_coarse'
    arg_types = ('virtual', 'state', 'iemaps', 'pbase' )
    geometry = [(Volume, 'virtual')]

    def __init__(self, name, sign, **kwargs):
        Term.__init__(self, name, sign,
                      function=terms.dw_mass_scalar_fine_coarse, **kwargs)
        
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
