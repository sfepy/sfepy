import numpy as nm

from sfepy.base.base import Struct
from sfepy.terms.terms import Term
from sfepy.terms.extmods import terms

class ContactInfo(Struct):
    """
    Various contact-related data of contact terms.
    """
    pass

class ContactTerm(Term):
    r"""
    Contact term with a penalty function.

    The penalty function is defined as :math:`\varepsilon_N \langle g_N(\ul{u})
    \rangle`, where :math:`\varepsilon_N` is the normal penalty parameter and
    :math:`\langle g_N(\ul{u}) \rangle` are the Macaulay's brackets of the gap
    function :math:`g_N(\ul{u})`.

    This term has a dynamic connectivity of DOFs in its region.

    :Definition:

    .. math::
        \int_{\Gamma_{c}} \varepsilon_N \langle g_N(\ul{u}) \rangle \ul{n}
        \ul{v}

    :Arguments:
        - material : :math:`\varepsilon_N`
        - virtual  : :math:`\ul{v}`
        - state    : :math:`\ul{u}`
    """
    name = 'dw_contact'
    arg_types = ('material', 'virtual', 'state')
    arg_shapes = {'material' : '.: 1',
                  'virtual' : ('D', 'state'), 'state' : 'D'}
    integration = 'surface'

    def __init__(self, *args, **kwargs):
        Term.__init__(self, *args, **kwargs)

        self.detect = 2
        self.ci = None

    def call_function(self, fargs):
        try:
            out, status = self.function(*fargs)

        except (RuntimeError, ValueError):
            terms.errclear()
            raise

        if status:
            terms.errclear()
            raise ValueError('term evaluation failed! (%s)' % self.name)

        return out, status

    def eval_real(self, shape, fargs, mode='eval', term_mode=None,
                  diff_var=None, **kwargs):
        out, status = self.call_function(fargs)
        if mode != 'weak':
            raise ValueError('unsupported evaluation mode! (%s)' % mode)

        return out, status

    @staticmethod
    def function(out_cc):
        return out_cc, 0

    def get_fargs(self, epss, virtual, state,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        geo, _ = self.get_mapping(virtual)

        region = self.region

        if self.ci is None:
            self.ci = ContactInfo()

        # Uses field connectivity (higher order nodes).
        sd = state.field.surface_data[region.name]

        ISN = state.field.efaces.T.copy()
        nsd = region.dim
        ngp = geo.n_qp
        neq = state.n_dof
        nsn = ISN.shape[0]

        fis = region.get_facet_indices()
        elementID = fis[:, 0].copy()
        segmentID = fis[:, 1].copy()

        n = len(elementID)
        IEN = state.field.econn

        # geo.bf corresponds to the field approximation.
        H = nm.asfortranarray(geo.bf[0, :, 0, :])

        gname = state.field.gel.name
        ib = 3 if gname == '3_8' else 0

        bqpkey = (self.integral.order, sd.bkey)
        bqp = state.field.qp_coors[bqpkey]
        gw = bqp.weights

        basis = state.field.create_basis_context()
        bfg3d = basis.evaluate(bqp.vals[ib], diff=True)
        bfg = bfg3d[:, :region.tdim - 1, state.field.efaces[ib]]
        if gname == '3_4':
            bfg = bfg[:, [1, 0], :]

        dH  = nm.asfortranarray(bfg.ravel().reshape(((nsd - 1) * ngp, nsn)))

        X = nm.asfortranarray(state.field.coors)
        Um = nm.asfortranarray(state().reshape((-1, nsd)))
        xx = nm.asfortranarray(X + Um)

        import sfepy.mechanics.extmods.ccontres as cc

        GPs = nm.empty((n*ngp, 2*nsd+6), dtype=nm.float64, order='F')

        longestEdge, GPs = cc.get_longest_edge_and_gps(GPs, neq,
                                                       elementID, segmentID,
                                                       ISN, IEN, H, xx)

        AABBmin, AABBmax = cc.get_AABB(xx, longestEdge, IEN, ISN,
                                       elementID, segmentID, neq);

        AABBmin = AABBmin - (0.5*longestEdge);
        AABBmax = AABBmax + (0.5*longestEdge);
        N = nm.ceil((AABBmax - AABBmin) / (0.5*longestEdge)).astype(nm.int32)

        head, next = cc.init_global_search(N, AABBmin, AABBmax, GPs[:,:nsd])

        npd = region.tdim - 1
        GPs = cc.evaluate_contact_constraints(GPs, ISN, IEN, N,
                                              AABBmin, AABBmax,
                                              head, next, xx,
                                              elementID, segmentID,
                                              npd, neq, longestEdge)

        Gc = nm.zeros(neq, dtype=nm.float64)

        activeGPs = GPs[:, 2*nsd+3]
        # gap = GPs[:, nsd + 2]
        # print activeGPs
        # print gap
        # print 'active:', activeGPs.sum()

        if diff_var is None:
            max_num = 1
            keyContactDetection = self.detect
            keyAssembleKc = 0

        else:
            max_num = 4 * (nsd * nsn)**2 * ngp * GPs.shape[0]
            keyContactDetection = self.detect
            keyAssembleKc = 1

        # print 'max_num:', max_num
        vals = nm.empty(max_num, dtype=nm.float64)
        rows = nm.empty(max_num, dtype=nm.int32)
        cols = nm.empty(max_num, dtype=nm.int32)

        aux = cc.assemble_contact_residual_and_stiffness(Gc, vals, rows, cols,
                                                         GPs, ISN, IEN,
                                                         X, Um, H, dH, gw,
                                                         activeGPs, neq, npd,
                                                         epss,
                                                         keyContactDetection,
                                                         keyAssembleKc)
        Gc, vals, rows, cols, num = aux
        # print 'true num:', num

        # Uncomment this to detect only in the 1. iteration.
        #self.detect = max(0, self.detect - 1)
        if diff_var is None:
            from sfepy.discrete.variables import create_adof_conn
            rows = nm.unique(create_adof_conn(nm.arange(len(Gc)), sd.econn,
                                              nsd, 0))
            out_cc = (Gc[rows], rows, state)

        else:
            out_cc = (vals[:num], rows[:num], cols[:num], state, state)

        return out_cc,
