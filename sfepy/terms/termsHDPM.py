from sfepy.terms.terms import *

class HDPMDiffusionVelocitySIntegratedTerm(Term):
    name = 'd_hdpm_surfdvel'
    arg_types = ('material', 'parameter')
    integration = 'surface_extra'

    function = staticmethod(terms.d_hdpm_surfdvel)
        
    def __call__( self, diff_var = None, chunk_size = None, **kwargs ):
        """? move surface pressure grad part into a cache ?"""
        mat, par = self.get_args( **kwargs )
        ap, sg = self.get_approximation(par)
        n_fa, n_qp = ap.get_s_data_shape( self.integral_name,
                                          self.region.name )[:2]
        shape = (chunk_size, 1, 1, 1)

        sd = ap.surface_data[self.region.name]

        vec = par()
        for out, chunk in self.char_fun( chunk_size, shape ):
            lchunk = self.char_fun.get_local_chunk()
            status = self.function( out, vec, 0,
                                    mat, sg, sd.fis, lchunk, ap.econn, chunk )
            out1 = nm.sum( out )
            yield out1, chunk, status
