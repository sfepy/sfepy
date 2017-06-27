import numpy as nm

from sfepy.base.base import Struct
from sfepy.terms.terms import Term

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

        self.ci = None

    @staticmethod
    def function(out, geo, fmode):
        out_qp = nm.zeros((out.shape[0], geo.n_qp) + out.shape[2:],
                          dtype=out.dtype)
        status = geo.integrate(out, nm.ascontiguousarray(out_qp))

        return status

    def get_fargs(self, epss, virtual, state,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        geo, _ = self.get_mapping(virtual)

        region = self.region

        if self.ci is None:
            self.ci = ContactInfo()

        print(region.name)
        print(region.shape)
        print(region.facets)
        print(region.get_facet_indices())
        print(state.field.surface_data[region.name].fis)

        print(geo)
        print(geo.normal)

        mesh_coors = self.region.domain.mesh.coors
        print(mesh_coors[region.vertices])

        # Uses field connectivity (higher order nodes).
        sd = state.field.surface_data[region.name]

        # Uses mesh connectivity.
        sdg = self.region.domain.surface_groups[region.name]

        print(mesh_coors[sdg.econn])

        qps = self.get_physical_qps()
        qp_coors = qps.values
        u_qp = self.get(state, 'val').reshape(qp_coors.shape)

        # Deformed QP coordinates.
        coors = u_qp + qp_coors

        print(coors)

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
        dH = nm.asfortranarray(bfg[0, ...])

        X = nm.asfortranarray(state.field.coors)
        Um = nm.asfortranarray(state().reshape((-1, nsd)))
        xx = nm.asfortranarray(X + Um)

        import sfepy.mechanics.extmods.ccontres as cc

        GPs = nm.empty((n*ngp, 2*nsd+6), dtype=nm.float64, order='F')

        longestEdge, GPs = cc.get_longest_edge_and_gps(GPs, neq,
                                                       elementID, segmentID,
                                                       ISN, IEN, H, xx)

        print longestEdge
        print X[IEN[elementID[0]]][ISN]

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

        print AABBmin, AABBmax, N
        head, next = cc.init_global_search(N, AABBmin, AABBmax, GPs[:,:nsd])

        print head, next

        npd = region.tdim - 1
        GPs = cc.evaluate_contact_constraints(GPs, ISN, IEN, N,
                                              AABBmin, AABBmax,
                                              head, next, xx,
                                              elementID, segmentID,
                                              npd, neq, longestEdge)
        print GPs[:, 4]

        #
        Gc = nm.zeros(neq, dtype=nm.float64)

        activeGPs = GPs[:, 2*nsd+3]
        keyContactDetection = 1
        keyAssembleKc = 1

        max_num = 4 * nsd * nsn * ngp * GPs.shape[0]
        print 'max_num:', max_num
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
        print Gc.mean(), num
        print GPs
        print 'true num:', num
        from sfepy.base.base import debug; debug()

        if diff_var is None:
            fmode = 0

        else:
            fmode = 1

        return geo, fmode
