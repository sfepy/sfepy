from sfepy.terms.terms import *

##
# 22.08.2006, c
class LinearTractionTerm( Term ):
    r"""
    :Description:
    Linear traction forces (weak form), where, depending on dimension of
    'material' argument, :math:`\ull{\sigma} \cdot \ul{n}` is
    :math:`\bar{p} \ull{I} \cdot \ul{n}` for a given scalar pressure,
    :math:`\ul{f}` for a traction vector, and itself for a stress tensor.

    :Definition:
    .. math::
        \int_{\Gamma} \ul{v} \cdot \ull{\sigma} \cdot \ul{n}

    :Arguments:
        material : :math:`\ull{\sigma}`,
        virtual  : :math:`\ul{v}`
    """
    name = 'dw_surface_ltr'
    arg_types = ('material', 'virtual')
    geometry = [(Surface, 'virtual')]
    dof_conn_type = 'surface'

    function = staticmethod(terms.dw_surface_ltr)

    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        """
        Should work in scalar, vector and tensor modes (tensor probably broken).
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

##        sg.str( sys.stdout, 0 )
##         print ap.bf[sd.face_type]
##         pause()

##         if ap.bf[sd.face_type].shape != gbf.shape:
##             raise NotImplementedError, 'tractions on P1 edges only!'
        for out, chunk in self.char_fun( chunk_size, shape ):
            lchunk = self.char_fun.get_local_chunk()
##             print out.shape, lchunk.shape
##             print traction.shape
            status = self.function( out, bf, gbf,
                                    traction, sg, lchunk )
##             print out
##             print nm.sum( out )
##             pause()
            yield out, lchunk, status

class SurfaceJumpTerm(Term):
    r"""
    :Description:
    Interface jump condition.
    
    :Definition:
    .. math::
        \int_{\Gamma} q (p_1 - p_2 - c)

    :Arguments:
        material : :math:`c`,
        virtual  : :math:`q`,
        state_1  : :math:`p_1`,
        state_2  : :math:`p_2`
    """
    name = 'dw_jump'
    arg_types = ('material', 'virtual', 'state_1', 'state_2')
    geometry = [(Surface, 'virtual'), (Surface, 'state_1'), (Surface, 'state_2')]
    dof_conn_type = 'surface'

    function = staticmethod(terms.dw_jump)

    def __call__(self, diff_var=None, chunk_size=None, **kwargs):
        coef, virtual, state1, state2 = self.get_args(**kwargs)
        ap, sg = virtual.get_approximation(self.get_current_group(), 'Surface')
        n_fa, n_qp, dim, n_fp = ap.get_s_data_shape(self.integral_name,
                                                    self.region.name)
        if diff_var is None:
            shape, mode = (chunk_size, 1, n_fp, 1), 0
        elif diff_var == self.get_arg_name('state_1'):
            shape, mode = (chunk_size, 1, n_fp, n_fp), 1
        elif diff_var == self.get_arg_name('state_2'):
            shape, mode = (chunk_size, 1, n_fp, n_fp), 2
        else:
            raise StopIteration

        sd = ap.surface_data[self.region.name]
        bf = ap.get_base(sd.face_type, 0, self.integral_name)
        
        ap1, sg1 = self.get_approximation(state1, kind='Surface')
        sd1 = ap1.surface_data[self.region.name]

        ap2, sg2 = self.get_approximation(state2, kind='Surface')
        sd2 = ap2.surface_data[self.region.name]

        for out, chunk in self.char_fun( chunk_size, shape ):
            lchunk = self.char_fun.get_local_chunk()
            status = self.function(out, coef, state1(), state2(),
                                   bf, sg, sd1.econn, sd2.econn, lchunk, mode)
##             print out
##             print nm.sum( out )
##             pause()
            yield out, lchunk, status
    
