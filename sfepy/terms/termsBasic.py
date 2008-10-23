from terms import *
from terms_base import VectorOrScalar, ScalarScalarTH
from utils import fix_mat_qp_shape

class IntegrateVolumeTerm( Term ):
    r""":definition: $\int_\Omega y$,  $\int_\Omega \ul{y}$"""
    name = 'di_volume_integrate'
    arg_types = ('parameter',)
    geometry = [(Volume, 'parameter')]
    use_caches = {'state_in_volume_qp' : [['parameter']]}

    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign )
        
    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        par, = self.get_args( **kwargs )
        ap, vg = par.get_approximation( self.get_current_group(), 'Volume' )
        n_el, n_qp, dim, n_ep = ap.get_v_data_shape( self.integral_name )

        field_dim = par.field.dim[0]
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
    r""":definition: $\int_\Omega q$"""

    name = 'dw_volume_integrate'
    arg_types = ('virtual',)
    geometry = [(Volume, 'virtual')]

    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign )
        
    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        virtual, = self.get_args( **kwargs )
        ap, vg = virtual.get_approximation( self.get_current_group(), 'Volume' )
        n_el, n_qp, dim, n_ep = ap.get_v_data_shape( self.integral_name )

        field_dim = par.field.dim[0]
        assert_( field_dim == 1 )
        
        if diff_var is None:
            shape = (chunk_size, 1, field_dim * n_ep, 1 )
            mode = 0
        else:
            raise StopIteration

        bf = ap.get_base( 'v', 0, self.integral_name )
        for out, chunk in self.char_fun( chunk_size, shape ):
            bf_t = nm.tile( bf.transpose( (0, 2, 1) ),
                            (chunk.shape[0], 1, 1, 1) )
            status = vg.integrate_chunk( out, bf_t, chunk )

            yield out, chunk, 0

## 24.04.2007, c
class IntegrateSurfaceTerm( Term ):
    r""":definition: $\int_\Gamma y$, for vectors: $\int_\Gamma \ul{y} \cdot \ul{n}$"""
    name = 'd_surface_integrate'
    arg_types = ('parameter',)
    geometry = [(Surface, 'parameter')]
    use_caches = {'state_in_surface_qp' : [['parameter']]}

    ##
    # 24.04.2007, c
    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign )
        self.dof_conn_type = 'surface'

    ##
    # c: 24.04.2007, r: 15.01.2008
    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        """
        Integrates over surface.
        """
        par, = self.get_args( **kwargs )
        ap, sg = par.get_approximation( self.get_current_group(), 'Surface' )
        shape = (chunk_size, 1, 1, 1)

        sd = ap.surface_data[self.region.name]

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
    r""":description: Volume $L^2(\Omega)$ dot product for both scalar and
    vector fields.
    :definition: $\int_\Omega p r$, $\int_\Omega \ul{u} \cdot \ul{w}$"""
    name = 'd_volume_dot'
    arg_types = ('parameter_1', 'parameter_2')
    geometry = [(Volume, 'parameter_1'), (Volume, 'parameter_2')]
    use_caches = {'state_in_volume_qp' : [['parameter_1'], ['parameter_2']]}

    ##
    # 26.09.2007, c
    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign )
        
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
    r""":description: Surface $L^2(\Gamma)$ dot product for both scalar and
    vector fields.
    :definition: $\int_\Gamma p r$, $\int_\Gamma \ul{u} \cdot \ul{w}$"""
    name = 'd_surface_dot'
    arg_types = ('parameter_1', 'parameter_2')
    geometry = [(Surface, 'parameter_1'), (Surface, 'parameter_2')]
    use_caches = {'state_in_surface_qp' : [['parameter_1'], ['parameter_2']]}

    ##
    # 09.10.2007, c
    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign )
        
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
    r""":definition: $\int_{\Gamma} q$"""

    name = 'dw_surface_integrate'
    arg_types = ('material', 'virtual',)
    geometry = [(Surface, 'virtual')]

    ##
    # 30.06.2008, c
    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign )

        self.dof_conn_type = 'surface'

    ##
    # 30.06.2008, c
    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        mat, virtual, = self.get_args( **kwargs )
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
            status = sg.integrate_chunk( out, bf_t, lchunk, 1 )

            out = out*mat
            yield out, lchunk, 0

##
# 16.07.2007, c
class VolumeTerm( Term ):
    r""":description: Volume of a domain. Uses approximation of the parameter
    variable.
    :definition: $\int_\Omega 1$"""
    name = 'd_volume'
    arg_types = ('parameter',)
    geometry = [(Volume, 'parameter')]
    use_caches = {'volume' : [['parameter']]}

    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign )
        
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

##
# c: 06.05.2008
class AverageVolumeMatTerm( Term ):
    r""":description: Material parameter $m$ averaged in elements. Uses
    approximation of $y$ variable.
    :definition: $\forall K \in \Tcal_h: \int_{T_K} m / \int_{T_K} 1$
    :arguments: material : $m$ (can have up to two dimensions),
    parameter : $y$, shape : shape of material parameter
    parameter, mode : 'const' or 'vertex' or 'element_avg'
    """
    name = 'de_volume_average_mat'
    arg_types = ('material', 'parameter', 'shape', 'mode')
    geometry = [(Volume, 'parameter')]
    use_caches = {'mat_in_qp' : [['material']]}

    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign )
        
    ##
    # c: 06.05.2008, r: 06.05.2008
    def prepare_data( self, chunk_size = None, **kwargs ):
        mat, par, mat_shape, mode = self.get_args( **kwargs )
        ap, vg = par.get_approximation( self.get_current_group(), 'Volume' )
        n_el, n_qp, dim, n_ep = ap.get_v_data_shape( self.integral_name )

        cache = self.get_cache( 'mat_in_qp', 0 )
        mat_qp = cache( 'matqp', self.get_current_group(), 0,
                       mat = mat, ap = ap,
                       assumed_shapes = [(1, n_qp) + mat_shape,
                                        (n_el, n_qp) + mat_shape],
                       mode_in = mode )

        mat_qp = fix_mat_qp_shape( mat_qp, chunk_size )
        shape = (chunk_size, 1) + mat_qp.shape[2:]

        return vg, mat_qp, shape

    ##
    # c: 06.05.2008, r: 14.07.2008
    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        vg, mat_qp, shape = self.prepare_data( chunk_size, **kwargs )

        for out, chunk in self.char_fun( chunk_size, shape ):
            status = vg.integrate_chunk( out, mat_qp[chunk], chunk )
            out1 = out / vg.variable( 2 )[chunk]
            yield out1, chunk, status

##
# c: 05.03.2008
class IntegrateVolumeMatTerm( AverageVolumeMatTerm ):
    r""":description: Integrate material parameter $m$ over a domain. Uses
    approximation of $y$ variable.
    :definition: $\int_\Omega m$
    :arguments: material : $m$ (can have up to two dimensions),
    parameter : $y$, shape : shape of material parameter
    parameter, mode : 'const' or 'vertex' or 'element_avg'
    """
    name = 'di_volume_integrate_mat'

    ##
    # c: 05.03.2008, r: 06.05.2008
    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        mat_shape, = self.get_args( ['shape'], **kwargs )
        vg, mat_qp, shape = self.prepare_data( chunk_size, **kwargs )

        for out, chunk in self.char_fun( chunk_size, shape ):
            status = vg.integrate_chunk( out, mat_qp[chunk], chunk )
            out1 = nm.sum( out, 0 )
            out1.shape = mat_shape
            yield out1, chunk, status

##
# c: 05.03.2008
class WDotProductVolumeTerm( VectorOrScalar, Term ):
    r""":description: Volume $L^2(\Omega)$ weighted dot product for both scalar
    and vector (not implemented in weak form!) fields. Can be evaluated. Can
    use derivatives.
    :definition: $\int_\Omega y q p$, $\int_\Omega y \ul{v} \cdot \ul{u}$,
    $\int_\Omega y p r$, $\int_\Omega y \ul{u} \cdot \ul{w}$
    :arguments: material : weight function $y$"""
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

        if mat.ndim == 1:
            mat = mat[...,nm.newaxis]

        n_el, n_qp, dim, n_ep = self.data_shape

        cache = self.get_cache( 'mat_in_qp', 0 )
        mat_qp = cache( 'matqp', self.get_current_group(), 0,
                       mat = mat, ap = ap,
                       assumed_shapes = [(1, n_qp, 1, 1),
                                        (n_el, n_qp, 1, 1)],
                       mode_in = None )

        mat_qp = fix_mat_qp_shape( mat_qp, chunk_size )

        bf = ap.get_base( 'v', 0, self.integral_name )

        return (vec_qp, bf, mat_qp, vg), shape, mode

    def dw_volume_wdot( self, out, vec_qp, bf, mat_qp, vg, chunk, mode ):
        if self.vdim > 1:
            raise NotImplementedError

        bf_t = bf.transpose( (0, 2, 1) )
        if mode == 0:
            vec = bf_t * mat_qp[chunk] * vec_qp[chunk]
        else:
            vec = bf_t * mat_qp[chunk] * bf
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

        if mat.ndim == 1:
            mat = mat[...,nm.newaxis]

        n_el, n_qp, dim, n_ep = self.data_shape

        cache = self.get_cache( 'mat_in_qp', 0 )
        mat_qp = cache( 'matqp', self.get_current_group(), 0,
                       mat = mat, ap = ap,
                       assumed_shapes = [(1, n_qp, 1, 1),
                                        (n_el, n_qp, 1, 1)],
                       mode_in = None )

        mat_qp = fix_mat_qp_shape( mat_qp, chunk_size )

        return (vec1_qp, vec2_qp, mat_qp, vg), (chunk_size, 1, 1, 1), 0

    def d_volume_wdot( self, out, vec1_qp, vec2_qp, mat_qp, vg, chunk ):

        if self.vdim > 1:
            vec = mat_qp[chunk] * nm.sum( vec1_qp[chunk] * vec2_qp[chunk],
                                          axis = -1 )
        else:
            vec = mat_qp[chunk] * vec1_qp[chunk] * vec2_qp[chunk]
        status = vg.integrate_chunk( out, vec, chunk )
        return status

    def set_arg_types( self ):
        if self.mode == 'weak':
            self.function = self.dw_volume_wdot
            use_method_with_name( self, self.get_fargs_weak, 'get_fargs' )
            self.use_caches = {'state_in_volume_qp' : [['state']],
                               'mat_in_qp' : [['material']]}
        else:
            self.function = self.d_volume_wdot
            use_method_with_name( self, self.get_fargs_eval, 'get_fargs' )
            self.use_caches = {'state_in_volume_qp' : [['parameter_1'],
                                                       ['parameter_2']],
                               'mat_in_qp' : [['material']]}

##
# c: 03.04.2008
class WDotSProductVolumeOperatorTHTerm( ScalarScalarTH, Term ):
    r""":description: Fading memory volume $L^2(\Omega)$ weighted dot product
    for scalar fields. Can use derivatives.
    :definition: $\int_\Omega \left [\int_0^t \Gcal(t-\tau) p(\tau)
    \difd{\tau} \right] q$"""
    name = 'dw_volume_wdot_scalar_th'
    arg_types = ('ts', 'material', 'virtual', 'state')
    geometry = [(Volume, 'virtual'), (Volume, 'state')]
    use_caches = {'state_in_volume_qp' : [['state', {'state' : (-1,-1)}]]}

    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign, terms.dw_volume_wdot_scalar )

    def get_fargs( self, diff_var = None, chunk_size = None, **kwargs ):
        ts, mats, virtual, state = self.get_args( **kwargs )
        ap, vg = virtual.get_approximation( self.get_current_group(), 'Volume' )

        self.set_data_shape( ap )
        shape, mode = self.get_shape( diff_var, chunk_size )

        if (ts.step == 0) and (mode == 0):
            raise StopIteration

        n_qp = self.data_shape[1]
        bf = ap.get_base( 'v', 0, self.integral_name )

        if mode == 1:
            aux = nm.array( [0], ndmin = 4, dtype = nm.float64 )
            mat_qp = mats[0][nm.newaxis,:,nm.newaxis].repeat( n_qp, 0 )
            return (ts.dt, aux, bf, mat_qp, vg), shape, mode

        else:
            cache = self.get_cache( 'state_in_volume_qp', 0 )
            def iter_kernel():
                for ii, mat in enumerate( mats ):
                    mat_qp = mat[nm.newaxis,:,nm.newaxis].repeat( n_qp, 0 )
                    vec_qp = cache( 'state', self.get_current_group(), ii,
                                    state = state, get_vector = self.get_vector )
                    yield ii, (ts.dt, vec_qp, bf, mat_qp, vg)
            return iter_kernel, shape, mode

class AverageVariableTerm( Term ):
    r""":description: Variable $y$ averaged in elements.
    :definition: vector of $\forall K \in \Tcal_h: \int_{T_K} y /
    \int_{T_K} 1$"""
    name = 'de_average_variable'
    arg_types = ('parameter',)
    geometry = [(Volume, 'parameter')]
    use_caches = {'state_in_volume_qp' : [['parameter']]}

    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign )
        
    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        par, = self.get_args( **kwargs )
        ap, vg = par.get_approximation( self.get_current_group(), 'Volume' )

        cache = self.get_cache( 'state_in_volume_qp', 0 )
        vec = cache( 'state', self.get_current_group(), 0,
                     state = par, get_vector = self.get_vector )
        vdim = vec.shape[-1]
        shape = (chunk_size, 1, vdim, 1)

        for out, chunk in self.char_fun( chunk_size, shape ):
            status = vg.integrate_chunk( out, vec[chunk], chunk )
            out1 = out / vg.variable( 2 )[chunk]
            yield out1, chunk, status
