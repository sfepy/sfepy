import numpy as nm

from sfepy.base.base import Struct
from sfepy.terms.terms import Term
from sfepy.terms.extmods import terms
import sfepy.mechanics.extmods.ccontres as cc

class ContactInfo(Struct):
    """
    Various contact-related data of contact terms.
    """
    def __init__(self, region, integral, geo, state):
        # Uses field connectivity (higher order nodes).
        sd = state.field.extra_data[f'sd_{region.name}']

        ISN = state.field.efaces.T.copy()
        nsd = region.dim
        ngp = geo.n_qp
        nsn = ISN.shape[0]

        fis = region.get_facet_indices()
        elementID = fis[:, 0].copy()
        segmentID = fis[:, 1].copy()

        # geo.bf corresponds to the field approximation.
        H = nm.asfortranarray(geo.bf[0, :, 0, :])

        gname = state.field.gel.name
        ib = 3 if gname == '3_8' else 0

        bqpkey = (integral.order, sd.bkey)
        state.field.create_bqp(region.name, integral)
        bqp = state.field.qp_coors[bqpkey]

        basis = state.field.create_basis_context()
        bfg3d = basis.evaluate(bqp.vals[ib], diff=True)
        bfg = bfg3d[:, :region.tdim - 1, state.field.efaces[ib]]
        if gname == '3_4':
            bfg = bfg[:, [1, 0], :]

        dH  = nm.asfortranarray(bfg.ravel().reshape(((nsd - 1) * ngp, nsn)))

        GPs = nm.empty((len(elementID)*ngp, 2*nsd+6),
                       dtype=nm.float64, order='F')

        Struct.__init__(self, name='contact_info_%s' % region.name,
                        sd=sd, nsd=nsd, ngp=ngp, neq=state.n_dof,
                        nsn=nsn, npd=region.tdim - 1,
                        elementID=elementID, segmentID=segmentID,
                        IEN=state.field.econn, ISN=ISN,
                        gw=bqp.weights, H=H, dH=dH, GPs=GPs)

    def update(self, xx):
        longestEdge, GPs = cc.get_longest_edge_and_gps(
            self.GPs, self.neq, self.elementID, self.segmentID,
            self.ISN, self.IEN, self.H, xx)

        AABBmin, AABBmax = cc.get_AABB(xx, longestEdge, self.IEN, self.ISN,
                                       self.elementID, self.segmentID, self.neq)

        AABBmin = AABBmin - (0.5*longestEdge);
        AABBmax = AABBmax + (0.5*longestEdge);
        N = nm.ceil((AABBmax - AABBmin) / (0.5*longestEdge)).astype(nm.int32)

        head, next = cc.init_global_search(N, AABBmin, AABBmax,
                                           GPs[:,:self.nsd])
        GPs = cc.evaluate_contact_constraints(
            GPs, self.ISN, self.IEN, N, AABBmin, AABBmax, head, next, xx,
            self.elementID, self.segmentID, self.npd, self.neq, longestEdge)

        return GPs

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
    integration = 'facet'

    def __init__(self, *args, **kwargs):
        Term.__init__(self, *args, **kwargs)

        self.detect = 2
        self.ci = None

    def get_contact_info(self, geo, state, init_gps=False):
        if self.ci is None:
            self.ci = ContactInfo(self.region, self.integral, geo, state)

        return self.ci

    def call_function(self, out, fargs):
        try:
            out, status = self.function(out, *fargs)

        except (RuntimeError, ValueError):
            terms.errclear()
            raise

        if status:
            terms.errclear()
            raise ValueError('term evaluation failed! (%s)' % self.name)

        return out, status

    def eval_real(self, shape, fargs, mode='eval', term_mode=None,
                  diff_var=None, **kwargs):
        if mode == 'weak':
            out, status = self.call_function(None, fargs)

        else:
            out = nm.empty(shape, dtype=nm.float64)
            status = self.call_function(out, fargs)

        return out, status

    @staticmethod
    def function_weak(out, out_cc):
        return out_cc, 0

    @staticmethod
    def integrate(out, val_qp, geo, fmode):
        if fmode == 2:
            out[:] = val_qp
            status = 0

        else:
            status = geo.integrate(out, val_qp, fmode)

        return out, status

    @staticmethod
    def function(out, fun, *args):
        return fun(out, *args)

    def get_fargs(self, epss, virtual, state,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        geo, _ = self.get_mapping(virtual)

        ci = self.get_contact_info(geo, state)

        X = nm.asfortranarray(state.field.coors)
        Um = nm.asfortranarray(state().reshape((-1, ci.nsd)))
        xx = nm.asfortranarray(X + Um)

        GPs = ci.update(xx)

        if mode == 'weak':
            Gc = nm.zeros(ci.neq, dtype=nm.float64)

            activeGPs = GPs[:, 2*ci.nsd+3]
            # gap = GPs[:, nsd + 2]
            # print activeGPs
            # print gap
            # print 'active:', activeGPs.sum()

            if diff_var is None:
                max_num = 1
                keyContactDetection = self.detect
                keyAssembleKc = 0

            else:
                max_num = 4 * (ci.nsd * ci.nsn)**2 * ci.ngp * GPs.shape[0]
                keyContactDetection = self.detect
                keyAssembleKc = 1

            # print 'max_num:', max_num
            vals = nm.empty(max_num, dtype=nm.float64)
            rows = nm.empty(max_num, dtype=nm.int32)
            cols = nm.empty(max_num, dtype=nm.int32)

            aux = cc.assemble_contact_residual_and_stiffness(
                Gc, vals, rows, cols, ci.GPs, ci.ISN, ci.IEN, X, Um,
                ci.H, ci.dH, ci.gw, activeGPs, ci.neq, ci.npd,
                epss, keyContactDetection, keyAssembleKc)
            Gc, vals, rows, cols, num = aux
            # print 'true num:', num

            # Uncomment this to detect only in the 1. iteration.
            #self.detect = max(0, self.detect - 1)
            if diff_var is None:
                from sfepy.discrete.variables import create_adof_conn
                rows = nm.unique(create_adof_conn(nm.arange(len(Gc)),
                                                  ci.sd.econn,
                                                  ci.nsd, 0))
                out_cc = (Gc[rows], rows, state)

            else:
                out_cc = (vals[:num], rows[:num], cols[:num], state, state)

            return self.function_weak, out_cc

        elif mode in ('el_avg', 'qp'):
            fmode = {'el_avg' : 1, 'qp' : 2}[mode]

            if term_mode == 'gap':
                gap = GPs[:, ci.nsd + 2].reshape(-1, ci.ngp, 1, 1)
                gap[gap > 0] = 0.0
                return self.integrate, gap, geo, fmode

            else:
                raise ValueError('unsupported term mode in %s! (%s)'
                                 % (self.name, term_mode))

        else:
            raise ValueError('unsupported evaluation mode in %s! (%s)'
                             % (self.name, mode))

    def get_eval_shape(self, epss, virtual, state,
                       mode=None, term_mode=None, diff_var=None, **kwargs):
        n_el, n_qp, dim, n_en, n_c = self.get_data_shape(state)

        if mode != 'qp':
            n_qp = 1

        return (n_el, n_qp, 1, 1), state.dtype
