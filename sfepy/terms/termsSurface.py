from sfepy.terms.terms import *

##
# c: 13.11.2007, r: 30.04.2008
def fix_traction_shape( mat, n_nod ):
    if nm.isscalar( mat ):
        mat = nm.tile( mat, (n_nod, 1) )
    else:
        mat = nm.array( mat, ndmin = 1 )
        if mat.ndim < 2:
            mat = mat[:,nm.newaxis]
        if mat.shape[0] != n_nod:
            if mat.shape[0] == 1:
                mat = nm.tile( mat, (n_nod, 1) )
            else:
                mat = nm.tile( mat[:,0], (n_nod, 1) )
    return nm.ascontiguousarray( mat )

##
# 22.08.2006, c
class LinearTractionTerm( Term ):
    r""":description: Linear traction forces (weak form), where,
    depending on dimension of 'material' argument, $\ull{\sigma} \cdot
    \ul{n}$ is $\bar{p} \ull{I} \cdot \ul{n}$ for a given scalar pressure,
    $\ul{f}$ for a traction vector, and itself for a stress tensor.
    :definition: $\int_{\Gamma} \ul{v} \cdot \ull{\sigma} \cdot \ul{n}$
    """
    name = 'dw_surface_ltr'
    arg_types = ('material', 'virtual')
    geometry = [(Surface, 'virtual')]

    # 11.10.2006
    def __init__( self, region, name = name, sign = 1 ):
        Term.__init__( self, region, name, sign, terms.dw_surface_ltr )
        self.dof_conn_type = 'surface'

    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        """
        Should work in scalar, vector and tensor modes (tensor probably broken).
        Tractions defined in vertices -> using 'vertex' subset of leconn
        """
        traction, virtual = self.get_args( **kwargs )
        ap, sg = virtual.get_approximation( self.get_current_group(), 'Surface' )
        n_fa, n_qp, dim, n_fp = ap.get_s_data_shape( self.integral_name,
                                                     self.region.name )
        if diff_var is None:
            shape = (chunk_size, 1, dim * n_fp, 1)
        else:
            raise StopIteration

        sd = ap.surface_data[self.region.name]
        bf = ap.get_base( sd.face_type, 0, self.integral_name )
        gbf = ap.get_base( sd.face_type, 0, self.integral_name,
                           from_geometry = True )

        traction = fix_traction_shape( traction, sd.nodes.shape[0] )
##        sg.str( sys.stdout, 0 )
##         print ap.bf[sd.face_type]
##         pause()

##         if ap.bf[sd.face_type].shape != gbf.shape:
##             raise NotImplementedError, 'tractions on P1 edges only!'
        n_ep = gbf.shape[2]
        leconn = sd.leconn[:,:n_ep].copy()
        for out, chunk in self.char_fun( chunk_size, shape ):
            lchunk = self.char_fun.get_local_chunk()
#            print out.shape, lchunk.shape
            status = self.function( out, bf, gbf,
                                    traction, sg, leconn, lchunk )
##             print out
##             print nm.sum( out )
##             pause()
            yield out, lchunk, status
