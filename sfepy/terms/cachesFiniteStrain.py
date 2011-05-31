from sfepy.terms.extmods import terms
from sfepy.terms.cache import DataCache

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

    def init_data(self, key, ckey, term, **kwargs):
        state, = self.get_args( **kwargs )

        n_el, n_qp, dim = state.get_data_shape(ckey[-1], term.integral)[:3]
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

    def update(self, key, term, ih, **kwargs):
        state, = self.get_args( **kwargs )
        ap, vg = term.get_approximation(state)

        ckey = self.get_key(term)

        try:
            self.function(self.data['F'][ckey][ih],
                          self.data['detF'][ckey][ih],
                          self.data['C'][ckey][ih],
                          self.data['trC'][ckey][ih],
                          self.data['in2C'][ckey][ih],
                          self.data['invC'][ckey][ih],
                          self.data['E'][ckey][ih],
                          state(), 0, vg, ap.econn)

        except RuntimeError:
            terms.errclear()
            raise ValueError

        self.valid['F'][ckey] = True
        self.valid['detF'][ckey] = True
        self.valid['C'][ckey] = True
        self.valid['trC'][ckey] = True
        self.valid['in2C'][ckey] = True
        self.valid['invC'][ckey] = True
        self.valid['E'][ckey] = True

class FiniteStrainSurfaceTLDataCache(DataCache):
    """
    Stores shared deformation-related data useful for the total Lagrangian
    formulation of finite strain elasticity for surface integrals.

    arguments:
      - state : displacements

    supported:
      - deformation gradient F
      - inverse deformation gradient F^{-1}
      - jacobian  J = det( F )

    All data are computed together, no matter which one is requested!
    """
    name = 'finite_strain_surface_tl'
    arg_types = ('state', 'data_shape')

    def __init__( self, name, arg_names, history_sizes = None ):

        keys = ['F', 'detF', 'invF']
        DataCache.__init__( self, name, arg_names, keys, history_sizes,
                            terms.dq_tl_finite_strain_surface )

    def init_data(self, key, ckey, term, **kwargs):
        state, data_shape = self.get_args( **kwargs )

        n_fa, n_qp, dim, n_ep = data_shape

        self.shapes = {
            'F' : (n_fa, n_qp, dim, dim),
            'invF' : (n_fa, n_qp, dim, dim),
            'detF' : (n_fa, n_qp, 1, 1),
        }

        DataCache.init_datas( self, ckey, self.shapes )

    def update(self, key, term, ih, **kwargs):
        state = self.get_args(**kwargs)[0]
        ap, sg = term.get_approximation(state)
        sd = ap.surface_data[term.region.name]

        ckey = self.get_key(term)

        try:
            self.function(self.data['F'][ckey][ih],
                          self.data['detF'][ckey][ih],
                          self.data['invF'][ckey][ih],
                          state(), 0, sg, sd.fis, ap.econn)

        except RuntimeError:
            terms.errclear()
            raise ValueError

        self.valid['F'][ckey] = True
        self.valid['detF'][ckey] = True
        self.valid['invF'][ckey] = True

class FiniteStrainULDataCache( DataCache ):
    """
    Stores shared deformation-related data useful for the updated Lagrangian
    formulation of finite strain elasticity.

    arguments:
      - state(s) : displacements u(t), u(t-1)

    supported:
      - relative deformation gradient F
      - jacobian  J = det( F )
      - left Cauchy-Green deformation tensor b = F F^T in vector (symmetric)
      storage
      - 1st invariant of b : tr( b ) 
      - 2nd invariant of b
     ???  - Green strain E = 1/2 (C - I)

    All data are computed together, no matter which one is requested!
    """
    name = 'finite_strain_ul'
    arg_types = ('state',)

    def __init__( self, name, arg_names, history_sizes = None ):

        keys = ['F', 'detF', 'B', 'trB', 'in2B', 'E']
        DataCache.__init__( self, name, arg_names, keys, history_sizes,
                            terms.dq_finite_strain_ul )

    def init_data(self, key, ckey, term, **kwargs):
        state, = self.get_args( **kwargs )

        n_el, n_qp, dim = state.get_data_shape(ckey[-1], term.integral)[:3]
        sym = dim * (dim + 1) / 2

        self.shapes = {
            'F' : (n_el, n_qp, dim, dim),
            'detF' : (n_el, n_qp, 1, 1),
            'B' : (n_el, n_qp, sym, 1),
            'trB' : (n_el, n_qp, 1, 1),
            'in2B' : (n_el, n_qp, 1, 1),
            'E' : (n_el, n_qp, sym, 1),
        }

        DataCache.init_datas( self, ckey, self.shapes )

    def update(self, key, term, ih, **kwargs):
        state, = self.get_args( **kwargs )
        ap, vg = term.get_approximation(state, geom0=True)

        ckey = self.get_key(term)

        try:
            self.function(self.data['F'][ckey][ih],
                          self.data['detF'][ckey][ih],
                          self.data['B'][ckey][ih],
                          self.data['trB'][ckey][ih],
                          self.data['in2B'][ckey][ih],
                          self.data['E'][ckey][ih],
                          state(), 0, vg, ap.econn)

        except RuntimeError:
            terms.errclear()
            raise ValueError

        self.valid['F'][ckey] = True
        self.valid['detF'][ckey] = True
        self.valid['B'][ckey] = True
        self.valid['trB'][ckey] = True
        self.valid['in2B'][ckey] = True
        self.valid['E'][ckey] = True
