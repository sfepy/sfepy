from sfepy.terms.terms import *

##
# 10.07.2007, c
class LinearPointSpringTerm( Term ):
    r"""
    :Description:
    Linear springs constraining movement of FE nodes in a region; to use as a
    relaxed Dirichlet boundary conditions.

    :Definition:
    .. math::
        \ul{f}^i = -k \ul{u}^i \quad \forall \mbox{ FE node } i \mbox{ in
        a region }

    :Arguments:
        material : :math:`k`,
        virtual  : :math:`\ul{v}`,
        state    : :math:`\ul{u}`
    """
    name = 'dw_point_lspring'
    arg_types = ('material', 'virtual', 'state')
    geometry = [(Point, 'virtual')]

    def __init__(self, name, sign, **kwargs):
        Term.__init__(self, name, sign, dof_conn_type='point', **kwargs)

    ##
    # 10.07.2007, c
    # 18.07.2007
    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        """TODO: projection to direction"""
        mat, virtual, state = self.get_args( **kwargs )
        ap, pg = virtual.get_approximation( self.get_current_group(), 'Point' )
        n_el, n_qp, dim, n_ep = ap.get_v_data_shape( self.integral_name )

        if self.char_fun.i_current > 0:
            raise StopIteration
            
        vec_u = state.get_state_in_region( self.region )
        n_nod = vec_u.shape[0]
##         print vec_u.shape
##         pause()

        stiffness = mat['stiffness']
        if diff_var is None:
            shape = (chunk_size, 1, dim, 1)
            for out, chunk in vector_chunk_generator( n_nod, chunk_size, shape ):
                out[:,0,:,0] = - stiffness * vec_u[chunk]
                yield out, chunk, 0

        elif diff_var == self.get_arg_name( 'state' ):
            shape = (chunk_size, 1, dim, dim)
            eye = nm.eye( dim, dim, dtype = nm.float64 )
            eye.shape = (1, 1) + eye.shape
            for out, chunk in vector_chunk_generator( n_nod, chunk_size, shape ):
                out[...] = - stiffness * eye
                yield out, chunk, 0

        else:
            raise StopIteration
