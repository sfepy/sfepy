import numpy as nm

from sfepy.base.base import output, Struct
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

class ContactIPCTerm(Term):
    r"""
    Contact term based on IPC Toolkit.

    This term has a dynamic connectivity of DOFs in its region.

    :Definition:

    .. math::
        \int_{\Gamma_{c}} \sum_{k \in C} \nabla b(d_k, \ul{u})

    :Arguments:
        - material_m : :math:`m`
        - material_k : :math:`k`
        - material_d : :math:`\hat{d}`
        - material_p : PSDProjectionMethod
        - virtual    : :math:`\ul{v}`
        - state      : :math:`\ul{u}`
    """
    name = 'dw_contact_ipc'
    arg_types = ('material_m', 'material_k', 'material_d', 'material_p',
                 'virtual', 'state')
    arg_shapes = {'material_m' : '.: 1', 'material_k' : '.: 1',
                  'material_d' : '.: 1', 'material_p' : '.: str',
                  'virtual' : ('D', 'state'), 'state' : 'D'}
    integration = 'facet'
    geometries = ['3_4', '3_8']

    def __init__(self, *args, **kwargs):
        Term.__init__(self, *args, **kwargs)

        import ipctk

        self.ipc = ipctk
        self.ci = None

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
    def function(out, fun, *args):
        return fun(out, *args)

    def get_contact_info(self, state):
        from sfepy.discrete.fem.mesh import Mesh
        from sfepy.mesh.mesh_tools import triangulate
        from sfepy.discrete.fem.geometry_element import create_geometry_elements
        from sfepy.discrete.variables import create_adof_conn

        if self.ci is not None:
            return self.ci

        tdim = self.region.kind_tdim
        if tdim not in (1, 2):
            raise ValueError(
                'contact region topological dimension must be 1 or 2!'
                f' (is {tdim})'
            )

        smesh = Mesh.from_region(self.region, self.region.domain.mesh,
                                 localize=True, is_surface=True)

        if '2_4' in smesh.descs:
            smesh = triangulate(smesh)

        gels = create_geometry_elements()

        cmesh = smesh.cmesh
        cmesh.set_local_entities(gels)
        cmesh.setup_entities()

        edges = cmesh.get_conn(1, 0).indices.astype(nm.int32)
        edges = edges.reshape((-1, 2))
        if tdim == 2:
            faces = cmesh.get_conn(2, 0).indices.astype(nm.int32)
            faces = faces.reshape((-1, 3))
            collision_mesh = self.ipc.CollisionMesh(smesh.coors, edges, faces)

        else:
            faces = None
            collision_mesh = self.ipc.CollisionMesh(smesh.coors, edges)

        nods = state.field.get_dofs_in_region(self.region, merge=True)

        econn = state.field.get_econn('facet', self.region)
        dofs = nm.unique(
            create_adof_conn(nm.arange(state.n_dof), econn, smesh.dim, 0)
        )

        self.ci = Struct(smesh=smesh, edges=edges, faces=faces, nods=nods,
                         econn=econn, dofs=dofs, dim=smesh.dim,
                         collision_mesh=collision_mesh,
                         prev_min_distance=None,
                         min_distance=None,
                         it=0,
                         ls_it=0,
                         barrier_potential=None,
                         barrier_stiffness=None,
                         max_barrier_stiffness=None,
                         e_grad_full=None,
                         bp_grad_norm=None,
                         e_grad_norm=None)
        return self.ci

    def get_fargs(self, avg_mass, stiffness, dhat, spd_projection,
                  virtual, state,
                  mode=None, term_mode=None, diff_var=None, ci=None,
                  **kwargs):
        if diff_var is not None:
            ci = self.ci

        if ci is None:
            output('WARNING: ci argument was not passed to dw_contact_ipc')
            ci = self.get_contact_info(state)

        collision_mesh = ci.collision_mesh

        uvec = state().reshape((-1, ci.dim))[ci.nods]

        vertices = collision_mesh.rest_positions + uvec
        collisions = self.ipc.Collisions()
        collisions.build(collision_mesh, vertices, dhat)

        ci.min_distance = collisions.compute_minimum_distance(
            collision_mesh, vertices,
        )
        ci.bbox_diagonal = self.ipc.world_bbox_diagonal_length(vertices)

        B = self.ipc.BarrierPotential(dhat)

        actual_stiffness = stiffness
        if actual_stiffness == 0.0: # Adaptive barrier stiffness
            if ci.it == 0:
                # This should be done at the start of each time step.
                bp_grad = B.gradient(collisions, collision_mesh, vertices)

                if ci.e_grad_full is not None:
                    e_grad = ci.e_grad_full[ci.dofs]

                else:
                    e_grad = nm.zeros_like(bp_grad)

                (ci.barrier_stiffness,
                 ci.max_barrier_stiffness) = self.ipc.initial_barrier_stiffness(
                     ci.bbox_diagonal,
                     B.barrier,
                     dhat, avg_mass,
                     e_grad, bp_grad,
                 )
                actual_stiffness = ci.barrier_stiffness
                ci.prev_min_distance = ci.min_distance
                ci.e_grad_norm = nm.linalg.norm(e_grad)
                ci.bp_grad_norm = nm.linalg.norm(bp_grad)

            elif (ci.ls_it == 0) and (diff_var is None):
                # Do not update during line-search!
                actual_stiffness = self.ipc.update_barrier_stiffness(
                    ci.prev_min_distance, ci.min_distance,
                    ci.max_barrier_stiffness, ci.barrier_stiffness,
                    ci.bbox_diagonal,
                )
                ci.barrier_stiffness = actual_stiffness
                ci.prev_min_distance = ci.min_distance

            else:
                actual_stiffness = ci.barrier_stiffness

        else:
            ci.barrier_stiffness = actual_stiffness
            ci.e_grad_norm = 0.0
            ci.bp_grad_norm = 0.0

        ci.barrier_potential = B(collisions, collision_mesh, vertices)

        if self.verbosity > 0:
            output('min. distance::', ci.min_distance)
            output('barrier potential:', ci.barrier_potential)
            output('barrier stiffness:', actual_stiffness)

        if diff_var is None:
            bp_grad = B.gradient(collisions, collision_mesh, vertices)

            out_cc = (actual_stiffness * bp_grad,
                      ci.dofs, state)

        else:
            proj = getattr(self.ipc.PSDProjectionMethod, spd_projection,
                           self.ipc.PSDProjectionMethod.NONE)
            bp_hess = B.hessian(
                collisions, collision_mesh, vertices,
                project_hessian_to_psd=proj,
            )

            ir, ic = bp_hess.nonzero()
            if len(ir):
                vals = bp_hess[ir, ic].A1

            else:
                vals = nm.array([], dtype=bp_hess.dtype)

            out_cc = (actual_stiffness * vals, ci.dofs[ir],
                      ci.dofs[ic], state, state)

        return self.function_weak, out_cc
