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

        # print(region.name)
        # print(region.shape)
        # print(region.facets)
        # print(region.get_facet_indices())
        # print(state.field.surface_data[region.name].fis)

        # print(geo)
        # print(geo.normal)

        # mesh_coors = self.region.domain.mesh.coors
        # print(mesh_coors[region.vertices])

        # Uses field connectivity (higher order nodes).
        sd = state.field.surface_data[region.name]

        # Uses mesh connectivity.
        # sdg = self.region.domain.surface_groups[region.name]

        # print(mesh_coors[sdg.econn])

        # qps = self.get_physical_qps()
        # qp_coors = qps.values
        # u_qp = self.get(state, 'val').reshape(qp_coors.shape)

        # Deformed QP coordinates.
        # coors = u_qp + qp_coors

        # print(coors)

        ISN = state.field.efaces.T.copy()
        nsd = region.dim
        ngp = geo.n_qp
        neq = state.n_dof
        nsn = ISN.shape[0]
        nes = ISN.shape[1]
        nen = state.field.gel.n_vertex

        fis = region.get_facet_indices()
        elementID = fis[:, 0].copy()
        segmentID = fis[:, 1].copy()

        n = len(elementID)
        IEN = state.field.econn
        # Need surface bf, bfg corresponding to field approximation here, not
        # geo...
        H = nm.asfortranarray(geo.bf[0, :, 0, :])
        ps = state.field.gel.surface_facet.poly_space
        gps, gw = self.integral.get_qp(state.field.gel.surface_facet.name)
        bfg = ps.eval_base(gps, diff=True)
        # ?? shape - try in 3D
        dH = nm.asfortranarray(bfg[:, 0 , :])

        X = nm.asfortranarray(state.field.coors)
        Um = nm.asfortranarray(state().reshape((-1, nsd)))
        xx = nm.asfortranarray(X + Um)

        import sfepy.mechanics.extmods.ccontres as cc

        GPs = nm.empty((n*ngp, 2*nsd+6), dtype=nm.float64, order='F')

        longestEdge, GPs = cc.get_longest_edge_and_gps(GPs, neq,
                                                       elementID, segmentID,
                                                       ISN, IEN, H, xx)

        # print longestEdge
        # print X[IEN[elementID[0]]][ISN]

        # gg = state.field.mappings[('Omega', 2, 'volume')]
        # Implement Region.get_diameters(dim).

        #print nm.sqrt(region.domain.get_element_diameters(state.field.region.cells, gg[0], 0)).max()
        #print (GPs[:, :2] == qp_coors).all()

        #GPs2 = nm.zeros((n*ngp, 2*nsd+6), dtype=nm.float64, order='F')
        #GPs2[:, :2] = qp_coors
        #GPs2[:, 2:4] = nm.repeat(fis, ngp, axis=0)
        #GPs2[:, 4] = nm.finfo(nm.float32).max

        #print (GPs == GPs2).all()
        #print nm.abs(GPs - GPs2).max()

        AABBmin, AABBmax = cc.get_AABB(xx, longestEdge, IEN, ISN,
                                       elementID, segmentID, neq);

        AABBmin = AABBmin - (0.5*longestEdge);
        AABBmax = AABBmax + (0.5*longestEdge);
        N = nm.ceil((AABBmax - AABBmin) / (0.5*longestEdge)).astype(nm.int32)
        N = nm.ones(nsd, dtype=nm.int32) # BUG workaround.

        # print AABBmin, AABBmax, N
        head, next = cc.init_global_search(N, AABBmin, AABBmax, GPs[:,:nsd])

        # print head, next

        npd = region.tdim - 1
        GPs = cc.evaluate_contact_constraints(GPs, ISN, IEN, N,
                                              AABBmin, AABBmax,
                                              head, next, xx,
                                              elementID, segmentID,
                                              npd, neq, longestEdge)
        # print GPs[:, 4]

        #
        Gc = nm.zeros(neq, dtype=nm.float64)

        activeGPs = GPs[:, 2*nsd+3]
        gap = GPs[:, nsd + 2]
        print activeGPs
        print gap
        print 'active:', activeGPs.sum()

        if diff_var is None:
            max_num = 1
            keyContactDetection = self.detect
            keyAssembleKc = 0

        else:
            max_num = 4 * (nsd * nsn)**2 * ngp * GPs.shape[0]
            keyContactDetection = self.detect
            keyAssembleKc = 1

        print 'max_num:', max_num
        vals = nm.empty(max_num, dtype=nm.float64)
        rows = nm.empty(max_num, dtype=nm.int32)
        cols = nm.empty(max_num, dtype=nm.int32)

        print '!!!detect', keyContactDetection
        aux = cc.assemble_contact_residual_and_stiffness(Gc, vals, rows, cols,
                                                         GPs, ISN, IEN,
                                                         X, Um, H, dH, gw,
                                                         activeGPs, neq, npd,
                                                         epss,
                                                         keyContactDetection,
                                                         keyAssembleKc)
        Gc, vals, rows, cols, num = aux
        # print Gc.mean()
        # print GPs
        print 'true num:', num
        print GPs.shape
        print activeGPs
        print gap

        # Uncomment this to detect only in the 1. iteration.
        #self.detect = max(0, self.detect - 1)
        if diff_var is None:
            from sfepy.discrete.variables import create_adof_conn
            rows = nm.unique(create_adof_conn(nm.arange(len(Gc)), sd.econn,
                                              nsd, 0))
            Gc = Gc[rows]

            eq = state.eq_map.eq
            erows = eq[rows]
            active = (erows >= 0)
            out_cc = (Gc[active], erows[active])

        else:
            vals, rows, cols = vals[:num], rows[:num], cols[:num]
            eq = state.eq_map.eq

            erows, ecols = eq[rows], eq[cols]
            active = (erows >= 0) & (ecols >= 0)
            out_cc = (vals[active], erows[active], ecols[active])

        return out_cc,
