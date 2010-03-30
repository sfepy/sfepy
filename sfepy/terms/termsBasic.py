from sfepy.terms.terms import *
from sfepy.terms.terms_base import VectorOrScalar, ScalarScalar, ScalarScalarTH

class IntegrateVolumeTerm( Term ):
    r"""
    :Description:
    Integrate a variable over a volume.
    

    :Definition:
    .. math::
        \int_\Omega y \mbox{ , } \int_\Omega \ul{y}
    """
    name = 'di_volume_integrate'
    arg_types = ('parameter',)
    geometry = [(Volume, 'parameter')]
    use_caches = {'state_in_volume_qp' : [['parameter']]}

    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        par, = self.get_args( **kwargs )
        ap, vg = par.get_approximation( self.get_current_group(), 'Volume' )
        n_el, n_qp, dim, n_ep = ap.get_v_data_shape( self.integral_name )

        field_dim = par.field.shape[0]
        shape = (chunk_size, 1, field_dim, 1)

        cache = self.get_cache( 'state_in_volume_qp', 0 )
        vec = cache( 'state', self.get_current_group(), 0,
                     state = par, get_vector = self.get_vector )

        for out, chunk in self.char_fun( chunk_size, shape ):
            status = vg.integrate_chunk( out, vec[chunk], chunk )
            out1 = nm.sum( out, 0 )
            out1.shape = (field_dim,)
            yield out1, chunk, status

class IntegrateVolumeOperatorTerm( Term ):
    r"""
    :Description:
    Volume integral of a test function.

    :Definition:
    .. math::
        \int_\Omega q
    """
    name = 'dw_volume_integrate'
    arg_types = ('virtual',)
    geometry = [(Volume, 'virtual')]

    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        virtual, = self.get_args( **kwargs )
        ap, vg = virtual.get_approximation( self.get_current_group(), 'Volume' )
        n_el, n_qp, dim, n_ep = ap.get_v_data_shape( self.integral_name )

        field_dim = virtual.field.shape[0]
        assert_( field_dim == 1 )
        
        if diff_var is None:
            shape = (chunk_size, 1, field_dim * n_ep, 1 )
        else:
            raise StopIteration

        bf = ap.get_base( 'v', 0, self.integral_name )
        for out, chunk in self.char_fun( chunk_size, shape ):
            bf_t = nm.tile( bf.transpose( (0, 2, 1) ),
                            (chunk.shape[0], 1, 1, 1) )
            bf_t = nm.ascontiguousarray(bf_t)
            status = vg.integrate_chunk( out, bf_t, chunk )
            yield out, chunk, 0

class IntegrateVolumeVariableOperatorTerm( Term ):
    r"""
    :Description:
    Volume integral of a test function with variable coefficient.


    :Definition:
    .. math::
        \int_\Omega c q
    """
    name = 'dw_volume_integrate_variable'
    arg_types = ('material', 'virtual',)
    geometry = [(Volume, 'virtual')]

    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        mat, virtual, = self.get_args( **kwargs )
        ap, vg = virtual.get_approximation( self.get_current_group(), 'Volume' )
        n_el, n_qp, dim, n_ep = ap.get_v_data_shape( self.integral_name )
        field_dim = virtual.field.shape[0]
        assert_( field_dim == 1 )
        
        if diff_var is None:
            shape = (chunk_size, 1, field_dim * n_ep, 1 )
        else:
            raise StopIteration
        
        bf = ap.get_base( 'v', 0, self.integral_name )
        for out, chunk in self.char_fun( chunk_size, shape ):
            bf_t = nm.tile( bf.transpose( (0, 2, 1) ),
                            (chunk.shape[0], 1, 1, 1) )
            bf_t = nm.ascontiguousarray(bf_t)
            val = mat[chunk] * bf_t
            status = vg.integrate_chunk( out, val, chunk )
            
            yield out, chunk, 0

## 24.04.2007, c
class IntegrateSurfaceTerm( Term ):
    r"""
    :Description:
    Integrate a variable over a surface.
    
    :Definition:
    .. math::
        \int_\Gamma y \mbox{ , for vectors: } \int_\Gamma \ul{y} \cdot
        \ul{n}
    """
    name = 'd_surface_integrate'
    arg_types = ('parameter',)
    geometry = [(Surface, 'parameter')]
    use_caches = {'state_in_surface_qp' : [['parameter']]}

    def __init__(self, name, sign, **kwargs):
        Term.__init__(self, name, sign, dof_conn_type='surface', **kwargs)

    ##
    # c: 24.04.2007, r: 15.01.2008
    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        """
        Integrates over surface.
        """
        par, = self.get_args( **kwargs )
        ap, sg = par.get_approximation( self.get_current_group(), 'Surface' )
        shape = (chunk_size, 1, 1, 1)

        cache = self.get_cache( 'state_in_surface_qp', 0 )
        vec = cache( 'state', self.get_current_group(), 0, state = par )
        
        for out, chunk in self.char_fun( chunk_size, shape ):
            lchunk = self.char_fun.get_local_chunk()
            status = sg.integrate_chunk( out, vec[lchunk], lchunk, 0 )

            out1 = nm.sum( out )
            yield out1, chunk, status

##
# 26.09.2007, c
class DotProductVolumeTerm( Term ):
    r"""
    :Description:
    Volume :math:`L^2(\Omega)` dot product for both scalar and vector
    fields.

    :Definition:
    .. math::
        \int_\Omega p r \mbox{ , } \int_\Omega \ul{u} \cdot \ul{w}
    """
    name = 'd_volume_dot'
    arg_types = ('parameter_1', 'parameter_2')
    geometry = [(Volume, 'parameter_1'), (Volume, 'parameter_2')]
    use_caches = {'state_in_volume_qp' : [['parameter_1'], ['parameter_2']]}

    ##
    # created:       26.09.2007
    # last revision: 13.12.2007
    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        par1, par2 = self.get_args( **kwargs )
        ap, vg = par1.get_approximation( self.get_current_group(), 'Volume' )
        shape = (chunk_size, 1, 1, 1)

        cache = self.get_cache( 'state_in_volume_qp', 0 )
        vec1 = cache( 'state', self.get_current_group(), 0,
                      state = par1, get_vector = self.get_vector )
        cache = self.get_cache( 'state_in_volume_qp', 1 )
        vec2 = cache( 'state', self.get_current_group(), 0,
                      state = par2, get_vector = self.get_vector )

        for out, chunk in self.char_fun( chunk_size, shape ):
            if vec1.shape[-1] > 1:
                vec = nm.sum( vec1[chunk] * vec2[chunk], axis = -1 )
            else:
                vec = vec1[chunk] * vec2[chunk]
            status = vg.integrate_chunk( out, vec, chunk )
            out1 = nm.sum( out )
            yield out1, chunk, status

##
# 09.10.2007, c
class DotProductSurfaceTerm( Term ):
    r"""
    :Description:
    Surface :math:`L^2(\Gamma)` dot product for both scalar and vector
    fields.

    :Definition:
    .. math::
        \int_\Gamma p r \mbox{ , } \int_\Gamma \ul{u} \cdot \ul{w}
    """
    name = 'd_surface_dot'
    arg_types = ('parameter_1', 'parameter_2')
    geometry = [(Surface, 'parameter_1'), (Surface, 'parameter_2')]
    use_caches = {'state_in_surface_qp' : [['parameter_1'], ['parameter_2']]}

    ##
    # c: 09.10.2007, r: 15.01.2008
    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        par1, par2 = self.get_args( **kwargs )
        ap, sg = par1.get_approximation( self.get_current_group(), 'Surface' )
        shape = (chunk_size, 1, 1, 1)

        cache = self.get_cache( 'state_in_surface_qp', 0 )
        vec1 = cache( 'state', self.get_current_group(), 0, state = par1 )
        cache = self.get_cache( 'state_in_surface_qp', 1 )
        vec2 = cache( 'state', self.get_current_group(), 0, state = par2 )

        for out, chunk in self.char_fun( chunk_size, shape ):
            lchunk = self.char_fun.get_local_chunk()
            if vec1.shape[-1] > 1:
                vec = nm.sum( vec1[lchunk] * vec2[lchunk], axis = -1 )
            else:
                vec = vec1[lchunk] * vec2[lchunk]

            status = sg.integrate_chunk( out, vec, lchunk, 0 )

            out1 = nm.sum( out )
            yield out1, chunk, status


##
# 30.06.2008, c
class IntegrateSurfaceOperatorTerm( Term ):
    r"""
    :Description:
    Surface integral of a test function.

    :Definition:
    .. math::
        \int_{\Gamma} q
    """
    name = 'dw_surface_integrate'
    arg_types = ('virtual',)
    geometry = [(Surface, 'virtual')]

    def __init__(self, name, sign, **kwargs):
        Term.__init__(self, name, sign, dof_conn_type='surface', **kwargs)

    ##
    # 30.06.2008, c
    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        virtual, = self.get_args( **kwargs )
        ap, sg = virtual.get_approximation( self.get_current_group(), 'Surface' )
        n_fa, n_qp, dim, n_fp = ap.get_s_data_shape( self.integral_name,
                                                     self.region.name )
        if diff_var is None:
            shape = (chunk_size, 1, n_fp, 1 )
        else:
            raise StopIteration

        sd = ap.surface_data[self.region.name]
        bf = ap.get_base( sd.face_type, 0, self.integral_name )

        for out, chunk in self.char_fun( chunk_size, shape ):
            lchunk = self.char_fun.get_local_chunk()
            bf_t = nm.tile( bf.transpose( (0, 2, 1) ), (chunk.shape[0], 1, 1, 1) )
            bf_t = nm.ascontiguousarray(bf_t)
            status = sg.integrate_chunk( out, bf_t, lchunk, 1 )

            yield out, lchunk, 0

class IntegrateSurfaceVariableOperatorTerm(Term):
    r"""
    :Description:
    Surface integral of a test function with variable coefficient.

    :Definition:
    .. math::
        \int_\Gamma c q
    """
    name = 'dw_surface_integrate_variable'
    arg_types = ('material', 'virtual')
    geometry = [(Surface, 'virtual')]

    def __init__(self, name, sign, **kwargs):
        Term.__init__(self, name, sign, dof_conn_type='surface', **kwargs)
        
    def __call__(self, diff_var=None, chunk_size=None, **kwargs):
        mat, virtual = self.get_args(**kwargs)
        ap, sg = virtual.get_approximation(self.get_current_group(), 'Surface')
        n_fa, n_qp, dim, n_fp = ap.get_s_data_shape(self.integral_name,
                                                    self.region.name)
        if diff_var is None:
            shape = (chunk_size, 1, n_fp, 1)
        else:
            raise StopIteration

        sd = ap.surface_data[self.region.name]
        bf = ap.get_base(sd.face_type, 0, self.integral_name)

        ac = nm.ascontiguousarray

        for out, chunk in self.char_fun( chunk_size, shape ):
            lchunk = self.char_fun.get_local_chunk()
            bf_t = nm.tile(bf.transpose((0, 2, 1) ), (chunk.shape[0], 1, 1, 1))
            bf_t = nm.ascontiguousarray(bf_t)
            val =  mat[lchunk] * bf_t
            if virtual.is_real:
                status = sg.integrate_chunk(out, val, lchunk, 1)
            else:
                status_r = sg.integrate_chunk(out, ac(val.real), lchunk, 1)
                out_imag = nm.zeros_like(out)
                status_i = sg.integrate_chunk(out_imag, ac(val.imag), lchunk, 1)

                status = status_r or status_i
                out = out + 1j * out_imag
                
            yield out, lchunk, status

##
# 16.07.2007, c
class VolumeTerm( Term ):
    r"""
    :Description:
    Volume of a domain. Uses approximation of the parameter variable.

    :Definition:
    .. math::
        \int_\Omega 1
    """
    name = 'd_volume'
    arg_types = ('parameter',)
    geometry = [(Volume, 'parameter')]
    use_caches = {'volume' : [['parameter']]}

    ##
    # created:       16.07.2007
    # last revision: 13.12.2007
    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        par, = self.get_args( **kwargs )
        shape = (1, 1, 1, 1)

        cache = self.get_cache( 'volume', 0 )
        volume = cache( 'volume', self.get_current_group(), 0,
                        region = self.char_fun.region, field = par.field )
        yield volume, 0, 0

class SurfaceTerm( Term ):
    r"""
    :Description:
    Surface of a domain. Uses approximation of the parameter variable.

    :Definition:
    .. math::
        \int_\Gamma 1
    """
    name = 'd_surface'
    arg_types = ('parameter',)
    geometry = [(Surface, 'parameter')]
    use_caches = {'surface' : [['parameter']]}

    def __init__(self, name, sign, **kwargs):
        Term.__init__(self, name, sign, dof_conn_type='surface', **kwargs)
        
    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        par, = self.get_args( **kwargs )
        shape = (1, 1, 1, 1)

        cache = self.get_cache( 'surface', 0 )
        surface = cache( 'surface', self.get_current_group(), 0,
                         region = self.char_fun.region, field = par.field )
        yield surface, 0, 0

class VolumeSurfaceTerm( Term ):
    r"""
    :Description:
    Volume of a domain, using a surface integral. Uses approximation of the
    parameter variable.

    :Definition:
    .. math::
        \int_\Gamma \ul{x} \cdot \ul{n}
    """
    name = 'd_volume_surface'
    arg_types = ('parameter',)
    geometry = [(Surface, 'parameter')]

    def __init__(self, name, sign, **kwargs):
        Term.__init__(self, name, sign, dof_conn_type='surface',
                      function=terms.d_volume_surface, **kwargs)

    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        par, = self.get_args( **kwargs )
        ap, sg = par.get_approximation( self.get_current_group(), 'Surface' )
        shape = (chunk_size, 1, 1, 1)

        sd = ap.surface_data[self.region.name]
        bf = ap.get_base( sd.face_type, 0, self.integral_name )
        coor = par.field.get_coor()
        for out, chunk in self.char_fun( chunk_size, shape ):
            lchunk = self.char_fun.get_local_chunk()
            status = self.function( out, coor, bf,
                                    sg, sd.econn.copy(), lchunk )

            out1 = nm.sum( out )
            yield out1, chunk, status

class SurfaceMomentTerm(Term):
    r"""
    :Description:
    Surface integral of the outer product of the unit outward normal
    :math:`\ul{n}` and the coordinate :math:`\ul{x}` shifted by :math:`\ul{x}_0`
    
    :Definition:
    .. math::
        \int_{\Gamma} \ul{n} (\ul{x} - \ul{x}_0)
    """
    name = 'di_surface_moment'
    arg_types = ('parameter', 'shift')
    geometry = [(Surface, 'parameter')]

    def __init__(self, name, sign, **kwargs):
        Term.__init__(self, name, sign, dof_conn_type='surface',
                      function=terms.di_surface_moment, **kwargs)

    def __call__(self, diff_var=None, chunk_size=None, **kwargs):
        par, shift = self.get_args(**kwargs)
        ap, sg = par.get_approximation(self.get_current_group(), 'Surface')
        n_fa, n_qp, dim, n_fp = ap.get_s_data_shape(self.integral_name,
                                                    self.region.name)
        shape = (chunk_size, 1, dim, dim)

        sd = ap.surface_data[self.region.name]
        bf = ap.get_base(sd.face_type, 0, self.integral_name)
        coor = par.field.get_coor() \
               - nm.asarray(shift, dtype=nm.float64)[None,:]
        for out, chunk in self.char_fun(chunk_size, shape):
            lchunk = self.char_fun.get_local_chunk()
            status = self.function(out, coor, bf,
                                   sg, sd.econn.copy(), lchunk)

            out1 = nm.sum(out, axis=0).squeeze()
            yield out1, chunk, status

##
# c: 06.05.2008
class AverageVolumeMatTerm( Term ):
    r"""
    :Description:
    Material parameter :math:`m` averaged in elements. Uses
    approximation of :math:`y` variable.

    :Definition:
    .. math::
        \forall K \in \Tcal_h: \int_{T_K} m / \int_{T_K} 1

    :Arguments:
    material : :math:`m` (can have up to two dimensions),
    parameter : :math:`y`,
    """
    name = 'de_volume_average_mat'
    arg_types = ('material', 'parameter')
    geometry = [(Volume, 'parameter')]

    ##
    # c: 06.05.2008, r: 06.05.2008
    def prepare_data( self, chunk_size = None, **kwargs ):
        mat, par = self.get_args( **kwargs )
        ap, vg = par.get_approximation( self.get_current_group(), 'Volume' )
        n_el, n_qp, dim, n_ep = ap.get_v_data_shape( self.integral_name )

        shape = (chunk_size, 1) + mat.shape[2:]

        return vg, mat, shape

    ##
    # c: 06.05.2008, r: 14.07.2008
    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        vg, mat, shape = self.prepare_data( chunk_size, **kwargs )

        for out, chunk in self.char_fun( chunk_size, shape ):
            status = vg.integrate_chunk( out, mat[chunk], chunk )
            out1 = out / vg.variable( 2 )[chunk]
            yield out1, chunk, status

##
# c: 05.03.2008
class IntegrateVolumeMatTerm( AverageVolumeMatTerm ):
    r"""
    :Description:
    Integrate material parameter :math:`m` over a domain. Uses
    approximation of :math:`y` variable.

    :Definition:
    .. math::
        \int_\Omega m

    :Arguments:
    material : :math:`m` (can have up to two dimensions),
    parameter : :math:`y`, shape : shape of material parameter
    """
    name = 'di_volume_integrate_mat'
    arg_types = ('material', 'parameter', 'shape')
    geometry = [(Volume, 'parameter')]

    def prepare_data( self, chunk_size = None, **kwargs ):
        mat, par, mat_shape = self.get_args( **kwargs )
        ap, vg = par.get_approximation( self.get_current_group(), 'Volume' )
        n_el, n_qp, dim, n_ep = ap.get_v_data_shape( self.integral_name )

        shape = (chunk_size, 1) + mat.shape[2:]

        return vg, mat, shape

    ##
    # c: 05.03.2008, r: 06.05.2008
    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        mat_shape, = self.get_args( ['shape'], **kwargs )
        vg, mat, shape = self.prepare_data( chunk_size, **kwargs )

        for out, chunk in self.char_fun( chunk_size, shape ):
            status = vg.integrate_chunk( out, mat[chunk], chunk )
            out1 = nm.sum( out, 0 )
            out1.shape = mat_shape
            yield out1, chunk, status

##
# c: 05.03.2008
class WDotProductVolumeTerm( VectorOrScalar, Term ):
    r"""
    :Description:
    Volume :math:`L^2(\Omega)` weighted dot product for both scalar
    and vector (not implemented in weak form!) fields. Can be evaluated. Can
    use derivatives.
    
    :Definition:
    .. math::
        \int_\Omega y q p \mbox{ , } \int_\Omega y \ul{v} \cdot \ul{u} \mbox{ , }
        \int_\Omega y p r \mbox{ , } \int_\Omega y \ul{u} \cdot \ul{w}

    :Arguments:
    material : weight function :math:`y`
    """
    name = 'dw_volume_wdot'
    arg_types = (('material', 'virtual', 'state'),
                 ('material', 'parameter_1', 'parameter_2'))
    geometry = ([(Volume, 'virtual')],
                [(Volume, 'parameter_1'), (Volume, 'parameter_2')])
    modes = ('weak', 'eval')

    def get_fargs_weak( self, diff_var = None, chunk_size = None, **kwargs ):
        mat, virtual, state = self.get_args( **kwargs )
        ap, vg = virtual.get_approximation( self.get_current_group(), 'Volume' )

        self.set_data_shape( ap )
        shape, mode = self.get_shape( diff_var, chunk_size )

        cache = self.get_cache( 'state_in_volume_qp', 0 )
        vec_qp = cache( 'state', self.get_current_group(), 0,
                        state = state, get_vector = self.get_vector )

        bf = ap.get_base( 'v', 0, self.integral_name )

        return (vec_qp, bf, mat, vg), shape, mode

    def dw_volume_wdot( self, out, vec_qp, bf, mat, vg, chunk, mode ):
        if self.vdim > 1:
            raise NotImplementedError

        bf_t = bf.transpose( (0, 2, 1) )
        if mode == 0:
            vec = bf_t * mat[chunk] * vec_qp[chunk]
        else:
            vec = bf_t * mat[chunk] * bf
        status = vg.integrate_chunk( out, vec, chunk )
        return status

    def get_fargs_eval( self, diff_var = None, chunk_size = None, **kwargs ):
        mat, par1, par2 = self.get_args( **kwargs )
        ap, vg = par1.get_approximation( self.get_current_group(), 'Volume' )

        self.set_data_shape( ap )

        cache = self.get_cache( 'state_in_volume_qp', 0 )
        vec1_qp = cache( 'state', self.get_current_group(), 0,
                         state = par1, get_vector = self.get_vector )
        cache = self.get_cache( 'state_in_volume_qp', 1 )
        vec2_qp = cache( 'state', self.get_current_group(), 0,
                         state = par2, get_vector = self.get_vector )

        return (vec1_qp, vec2_qp, mat, vg), (chunk_size, 1, 1, 1), 0

    def d_volume_wdot( self, out, vec1_qp, vec2_qp, mat, vg, chunk ):

        if self.vdim > 1:
            vec = mat[chunk] * nm.sum( vec1_qp[chunk] * vec2_qp[chunk],
                                          axis = -1 )
        else:
            vec = mat[chunk] * vec1_qp[chunk] * vec2_qp[chunk]
        status = vg.integrate_chunk( out, vec, chunk )
        return status

    def set_arg_types( self ):
        if self.mode == 'weak':
            self.function = self.dw_volume_wdot
            use_method_with_name( self, self.get_fargs_weak, 'get_fargs' )
            self.use_caches = {'state_in_volume_qp' : [['state']]}
        else:
            self.function = self.d_volume_wdot
            use_method_with_name( self, self.get_fargs_eval, 'get_fargs' )
            self.use_caches = {'state_in_volume_qp' : [['parameter_1'],
                                                       ['parameter_2']]}

##
# c: 03.04.2008
class WDotSProductVolumeOperatorTHTerm( ScalarScalarTH, Term ):
    r"""
    :Description:
    Fading memory volume :math:`L^2(\Omega)` weighted dot product for
    scalar fields. Can use derivatives.

    :Definition:
    .. math::
        \int_\Omega \left [\int_0^t \Gcal(t-\tau) p(\tau) \difd{\tau} \right] q
    """
    name = 'dw_volume_wdot_scalar_th'
    arg_types = ('ts', 'material', 'virtual', 'state')
    geometry = [(Volume, 'virtual'), (Volume, 'state')]
    use_caches = {'state_in_volume_qp' : [['state', {'state' : (-1,-1)}]]}

    function = staticmethod(terms.dw_volume_wdot_scalar)

    def get_fargs( self, diff_var = None, chunk_size = None, **kwargs ):
        ts, mats, virtual, state = self.get_args( **kwargs )
        ap, vg = virtual.get_approximation( self.get_current_group(), 'Volume' )

        self.set_data_shape( ap )
        shape, mode = self.get_shape( diff_var, chunk_size )

        if (ts.step == 0) and (mode == 0):
            raise StopIteration

        n_el, n_qp = self.data_shape[:2]
        bf = ap.get_base( 'v', 0, self.integral_name )

        if mode == 1:
            aux = nm.array( [0], ndmin = 4, dtype = nm.float64 )
            mat = mats[0]
            mat = nm.tile(mat, (n_el, n_qp, 1, 1))
            return (ts.dt, aux, bf, mat, vg), shape, mode

        else:
            cache = self.get_cache( 'state_in_volume_qp', 0 )
            def iter_kernel():
                for ii, mat in enumerate( mats ):
                    vec_qp = cache( 'state', self.get_current_group(), ii,
                                    state = state, get_vector = self.get_vector )
                    mat = nm.tile(mat, (n_el, n_qp, 1, 1))
                    yield ii, (ts.dt, vec_qp, bf, mat, vg)
            return iter_kernel, shape, mode

class WDotSProductVolumeOperatorETHTerm( ScalarScalar, Term ):
    r"""
    :Description:
    Fading memory volume :math:`L^2(\Omega)` weighted dot product for
    scalar fields. This term has the same definition as
    dw_volume_wdot_scalar_th, but assumes an exponential approximation of
    the convolution kernel resulting in much higher efficiency. Can use
    derivatives.

    :Definition:
    .. math::
        \int_\Omega \left [\int_0^t \Gcal(t-\tau) p(\tau) \difd{\tau} \right] q

    :Arguments:
    material_0 : :math:`\Gcal(0)`,
    material_1 : :math:`\exp(-\lambda \Delta t)` (decay at :math:`t_1`)
    """
    name = 'dw_volume_wdot_scalar_eth'
    arg_types = ('ts', 'material_0', 'material_1', 'virtual', 'state')
    geometry = [(Volume, 'virtual'), (Volume, 'state')]
    use_caches = {'state_in_volume_qp' : [['state']],
                  'exp_history' : [['material_0', 'material_1', 'state']]}

    function = staticmethod(terms.dw_volume_wdot_scalar)

    def get_fargs( self, diff_var = None, chunk_size = None, **kwargs ):
        ts, mat0, mat1, virtual, state = self.get_args( **kwargs )
        ap, vg = virtual.get_approximation( self.get_current_group(), 'Volume' )

        self.set_data_shape( ap )
        shape, mode = self.get_shape( diff_var, chunk_size )

        bf = ap.get_base( 'v', 0, self.integral_name )
        if diff_var is None:
            cache = self.get_cache( 'state_in_volume_qp', 0 )
            vec_qp = cache( 'state', self.get_current_group(), 0,
                            state = state, get_vector = self.get_vector )

            cache = self.get_cache('exp_history', 0)
            increment = cache('increment', self.get_current_group(), 0,
                              decay=mat1, values=vec_qp)
            history = cache('history', self.get_current_group(), 0)

            fargs = (ts.dt, history + increment, bf, mat0, vg)
            if ts.step == 0: # Just init the history in step 0.
                raise StopIteration

        else:
            aux = nm.array([0], ndmin=4, dtype=nm.float64)
            fargs = (ts.dt, aux, bf, mat0, vg)

        return fargs, shape, mode

class AverageVariableTerm( Term ):
    r"""
    :Description:
    Variable :math:`y` averaged in elements.

    :Definition:
    .. math::
        \mbox{vector of } \forall K \in \Tcal_h: \int_{T_K} y /
        \int_{T_K} 1
    """
    name = 'de_average_variable'
    arg_types = ('parameter',)
    geometry = [(Volume, 'parameter')]
    use_caches = {'state_in_volume_qp' : [['parameter']]}

    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        par, = self.get_args( **kwargs )
        ap, vg = par.get_approximation( self.get_current_group(), 'Volume' )

        cache = self.get_cache( 'state_in_volume_qp', 0 )
        vec = cache( 'state', self.get_current_group(), 0,
                     state = par, get_vector = self.get_vector )
        vdim = vec.shape[2]
        shape = (chunk_size, 1, vdim, 1)

        for out, chunk in self.char_fun( chunk_size, shape ):
            status = vg.integrate_chunk( out, vec[chunk], chunk )
            out1 = out / vg.variable( 2 )[chunk]
            yield out1, chunk, status

class StateVQTerm(Term):
    r"""
    :Description:
    State interpolated into volume quadrature points.

    :Definition:
    .. math::
        \ul{u}|_{qp} \mbox{ , } p|_{qp}
    """
    name = 'dq_state_in_volume_qp'
    arg_types = ('state',)
    geometry = [(Volume, 'state')]

    function = staticmethod(terms.dq_state_in_qp)

    def __call__(self, diff_var=None, chunk_size=None, **kwargs):
        """Ignores chunk_size."""
        state, = self.get_args(**kwargs)
        ap, vg = state.get_approximation(self.get_current_group(), 'Volume')
        n_el, n_qp, dim, n_ep = ap.get_v_data_shape(self.integral_name)

        if diff_var is None:
            shape = (n_el, n_qp, state.dpn, 1)
        else:
            raise StopIteration

        vec = self.get_vector(state)
        bf = ap.get_base('v', 0, self.integral_name)

        out = nm.empty(shape, dtype=nm.float64)
        self.function(out, vec, 0, bf, ap.econn)

        yield out, nm.arange(n_el, dtype=nm.int32), 0

class StateSQTerm(Term):
    r"""
    :Description:
    State interpolated into surface quadrature points.

    :Definition:
    .. math::
        \ul{u}|_{qp} \mbox{ , } p|_{qp}
    """
    name = 'dq_state_in_surface_qp'
    arg_types = ('state',)
    geometry = [(Surface, 'state')]

    function = staticmethod(terms.dq_state_in_qp)

    def __call__(self, diff_var=None, chunk_size=None, **kwargs):
        """Ignores chunk_size."""
        state, = self.get_args(**kwargs)
        ap, sg = virtual.get_approximation(self.get_current_group(), 'Surface')
        n_fa, n_qp, dim, n_fp = ap.get_s_data_shape(self.integral_name,
                                                    self.region.name)

        if diff_var is None:
            shape = (chunk_size, n_qp, state.dpn, 1)
        else:
            raise StopIteration

        vec = self.get_vector(state)

        sd = ap.surface_data[self.region.name]
        bf = ap.get_base(sd.face_type, 0, self.integral_name)

        out = nm.empty(shape, dtype=nm.float64)
        self.function(out, vec, 0, bf, sd.econn)

        yield out, nm.arange(n_fa, dtype=nm.int32), 0
