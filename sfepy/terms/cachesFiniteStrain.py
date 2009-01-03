import numpy as nm
from sfepy.terms.extmods import terms
from sfepy.terms.cache import DataCache
from sfepy.base.base import pause, debug

class FiniteStrainTLDataCache( DataCache ):
    """
    Stores shared deformation-related data useful for the total Lagrangian
    formulation of finite strain elasticity.

    arguments:
      - state : displacements

    supported:
      - deformation gradient F
      - jacobian  J = det( F )
      - right Cauchy-Green deformation tensor C = F^T F in vector (symmetric)
      storage
      - 1st invariant of C : tr( C ) 
      - 2nd invariant of C
      - C^{-1} in vector (symmetric) storage
      - Green strain E = 1/2 (C - I)

    All data are computed together, no matter which one is requested!
    """
    name = 'finite_strain_tl'
    arg_types = ('state',)

    def __init__( self, name, arg_names, history_sizes = None ):

        keys = ['F', 'detF', 'C', 'trC', 'in2C', 'invC', 'E']
        DataCache.__init__( self, name, arg_names, keys, history_sizes,
                            terms.dq_finite_strain_tl )

    def init_data( self, key, ckey, **kwargs ):
        state, = self.get_args( **kwargs )

        n_el, n_qp, dim = state.get_data_shapes( ckey )[:3]
        sym = dim * (dim + 1) / 2

        self.shapes = {
            'F' : (n_el, n_qp, dim, dim),
            'detF' : (n_el, n_qp, 1, 1),
            'C' : (n_el, n_qp, sym, 1),
            'trC' : (n_el, n_qp, 1, 1),
            'in2C' : (n_el, n_qp, 1, 1),
            'invC' : (n_el, n_qp, sym, 1),
            'E' : (n_el, n_qp, sym, 1),
        }

        DataCache.init_datas( self, ckey, self.shapes )

    def update( self, key, group_indx, ih, **kwargs ):
##         print 'update!', key, group_indx, ih
        state, = self.get_args( **kwargs )
        ap, vg = state.get_approximation( group_indx, 'Volume' )

        ckey = self.g_to_c( group_indx )

        self.function( self.data['F'][ckey][ih],
                       self.data['detF'][ckey][ih],
                       self.data['C'][ckey][ih],
                       self.data['trC'][ckey][ih],
                       self.data['in2C'][ckey][ih],
                       self.data['invC'][ckey][ih],
                       self.data['E'][ckey][ih],
                       state(), 0, vg, ap.econn )
        
        self.valid['F'][ckey] = True
        self.valid['detF'][ckey] = True
        self.valid['C'][ckey] = True
        self.valid['trC'][ckey] = True
        self.valid['in2C'][ckey] = True
        self.valid['invC'][ckey] = True
        self.valid['E'][ckey] = True
