import numpy as nm

from sfepy.base.base import Struct, assert_
from sfepy.discrete.fem.mappings import VolumeMapping, SurfaceMapping
from poly_spaces import PolySpace
from fe_surface import FESurface

def set_mesh_coors(domain, fields, coors, update_fields=False, actual=False,
                   clear_all=True):
    if actual:
        domain.mesh.coors_act = coors.copy()
    else:
        domain.mesh.coors = coors.copy()
        domain.cmesh.coors[:] = coors

    if update_fields:
        for field in fields.itervalues():
            field.setup_coors(coors)
            field.clear_mappings(clear_all=clear_all)

def eval_nodal_coors(coors, mesh_coors, region, poly_space, geom_poly_space,
                     econn, ig, only_extra=True):
    """
    Compute coordinates of nodes corresponding to `poly_space`, given
    mesh coordinates and `geom_poly_space`.
    """
    if only_extra:
        iex = (poly_space.nts[:,0] > 0).nonzero()[0]
        if iex.shape[0] == 0: return

        qp_coors = poly_space.node_coors[iex, :]
        econn = econn[:, iex].copy()

    else:
        qp_coors = poly_space.node_coors

    ##
    # Evaluate geometry interpolation base functions in (extra) nodes.
    bf = geom_poly_space.eval_base(qp_coors)
    bf = bf[:,0,:].copy()

    ##
    # Evaluate extra coordinates with 'bf'.
    group = region.domain.groups[ig]
    cells = region.get_cells(ig)

    ecoors = nm.dot(bf, mesh_coors[group.conn[cells]])
    coors[econn] = nm.swapaxes(ecoors, 0, 1)


##
# 04.08.2005, c
def _interp_to_faces( vertex_vals, bfs, faces ):
    dim = vertex_vals.shape[1]
    n_face = faces.shape[0]
    n_qp = bfs.shape[0]
    
    faces_vals = nm.zeros( (n_face, n_qp, dim), nm.float64 )
    for ii, face in enumerate( faces ):
        vals = vertex_vals[face,:dim]
        faces_vals[ii,:,:] = nm.dot( bfs[:,0,:], vals )

    return( faces_vals )

class Interpolant( Struct ):
    """A simple wrapper around PolySpace."""

    def __init__(self, name, gel, space='H1', base='lagrange',
                 approx_order=1, force_bubble=False):
        self.name = name
        self.gel = gel

        self.poly_spaces = poly_spaces = {}
        poly_spaces['v'] = PolySpace.any_from_args(name, gel, approx_order,
                                                   base=base,
                                                   force_bubble=force_bubble)
        gel = gel.surface_facet
        if gel is not None:
            ps = PolySpace.any_from_args(name, gel, approx_order,
                                         base=base,
                                         force_bubble=False)
            skey = 's%d' % ps.n_nod
            poly_spaces[skey] = ps

    def describe_nodes( self ):
        ps = self.poly_spaces['v']
        node_desc = ps.describe_nodes()

        return node_desc

    ##
    # 16.11.2007, c
    def get_n_nodes( self ):
        nn = {}
        for key, ps in self.poly_spaces.iteritems():
            nn[key] = ps.nodes.shape[0]
        return nn

    def get_geom_poly_space(self, key):
        if key == 'v':
            ps = self.gel.interp.poly_spaces['v']

        elif key[0] == 's':
            n_v = self.gel.surface_facet.n_vertex
            ps = self.gel.interp.poly_spaces['s%d' % n_v]

        else:
            raise ValueError('bad polynomial space key! (%s)' % key)

        return ps

class SurfaceInterpolant(Interpolant):
    """
    Like Interpolant, but for use with SurfaceField and
    SurfaceApproximation.
    """

    def __init__(self, name, gel, space='H1', base='lagrange',
                 approx_order=1, force_bubble=False):
        Interpolant.__init__(self, name, gel, space=space, base=base,
                             approx_order=approx_order,
                             force_bubble=force_bubble)

        # Make alias 'v' <-> 's#'.
        ps = self.poly_spaces['v']
        self.poly_spaces['s%d' % ps.n_nod] = ps

    def get_geom_poly_space(self, key):
        assert_(key[0] == 's')

        ps = self.gel.interp.poly_spaces['v']

        return ps

##
# 18.07.2006, c
class Approximation( Struct ):
    ##
    # 18.07.2006, c
    # 10.10.2006
    # 11.07.2007
    # 17.07.2007
    def __init__(self, name, interp, region, ig, is_surface=False):
        """interp, region are borrowed."""

        self.name = name
        self.interp = interp
        self.region = region
        self.ig = ig
        self.is_surface = is_surface
        self.surface_data = {}
        self.edge_data = {}
        self.point_data = {}
        self.n_ep = self.interp.get_n_nodes()
        self.ori = None

        self.clear_qp_base()

    def eval_extra_coor(self, coors, mesh_coors):
        """
        Compute coordinates of extra nodes.
        """
        gps = self.interp.gel.interp.poly_spaces['v']
        ps = self.interp.poly_spaces['v']

        eval_nodal_coors(coors, mesh_coors, self.region, ps, gps, self.econn, self.ig)

    ##
    # c: 05.09.2006, r: 09.05.2008
    def setup_surface_data( self, region ):
        """nodes[leconn] == econn"""
        """nodes are sorted by node number -> same order as region.vertices"""
        sd = FESurface('surface_data_%s' % region.name, region,
                       self.efaces, self.econn, self.ig)
        self.surface_data[region.name] = sd
        return sd

    ##
    # 11.07.2007, c
    def setup_point_data( self, field, region ):
        conn = field.get_dofs_in_region(region, merge=True, igs=region.igs)
##         conn = [nods]\
##                + [nm.empty( (0,), dtype = nm.int32 )]\
##                * (len( region.igs ) - 1)

        conn.shape += (1,)
        self.point_data[region.name] = conn

    def get_connectivity(self, region, integration, is_trace=False):
        """
        Return the DOF connectivity for the given geometry type.

        Parameters
        ----------
        region : Region instance
            The region, used to index surface and volume connectivities.
        integration : one of ('volume', 'plate', 'surface', 'surface_extra')
            The term integration type.
        """
        if integration == 'surface':
            sd = self.surface_data[region.name]
            conn = sd.get_connectivity(self.is_surface, is_trace=is_trace)

        elif integration in ('volume', 'plate', 'surface_extra'):
            if region.name == self.region.name:
                conn = self.econn

            else:
                aux = integration in ('volume', 'plate')
                cells = region.get_cells(self.ig, true_cells_only=aux)
                conn = nm.take(self.econn, cells.astype(nm.int32), axis=0)

        else:
            raise ValueError('unsupported term integration! (%s)' % integration)

        return conn

    def get_poly_space(self, key, from_geometry=False):
        """
        Get the polynomial space.

        Parameters
        ----------
        key : 'v' or 's?'
            The key denoting volume or surface.
        from_geometry : bool
            If True, return the polynomial space for affine geometrical
            interpolation.

        Returns
        -------
        ps : PolySpace instance
            The polynomial space.
        """
        if from_geometry:
            ps = self.interp.get_geom_poly_space(key)

        else:
            ps = self.interp.poly_spaces[key]

        return ps

    def clear_qp_base(self):
        """
        Remove cached quadrature points and base functions.
        """
        self.qp_coors = {}
        self.bf = {}

    def get_qp(self, key, integral):
        """
        Get quadrature points and weights corresponding to the given key
        and integral. The key is 'v' or 's#', where # is the number of
        face vertices.
        """
        qpkey = (integral.name, key)

        if not self.qp_coors.has_key(qpkey):
            interp = self.interp
            if (key[0] == 's'):
                dim = interp.gel.dim - 1
                n_fp = interp.gel.surface_facet.n_vertex
                geometry = '%d_%d' % (dim, n_fp)

            else:
                geometry = interp.gel.name

            vals, weights = integral.get_qp(geometry)
            self.qp_coors[qpkey] = Struct(vals=vals, weights=weights)

        return self.qp_coors[qpkey]

    def get_base(self, key, derivative, integral, iels=None,
                 from_geometry=False, base_only=True):
        qp = self.get_qp(key, integral)

        ps = self.get_poly_space(key, from_geometry=from_geometry)

        _key = key if not from_geometry else 'g' + key
        bf_key = (integral.name, _key, derivative)

        if not self.bf.has_key(bf_key):
            if (iels is not None) and (self.ori is not None):
                ori = self.ori[iels]

            else:
                ori = self.ori

            self.bf[bf_key] = ps.eval_base(qp.vals, diff=derivative, ori=ori)

        if base_only:
            return self.bf[bf_key]
        else:
            return self.bf[bf_key], qp.weights

    def describe_geometry(self, field, gtype, region, integral,
                          return_mapping=False):
        """
        Compute jacobians, element volumes and base function derivatives
        for Volume-type geometries (volume mappings), and jacobians,
        normals and base function derivatives for Surface-type
        geometries (surface mappings).

        Notes
        -----
        - volume mappings can be defined on a part of an element group,
          although the field has to be defined always on the whole group.
        - surface mappings are defined on the surface region
        - surface mappings require field order to be > 0
        """
        domain = field.domain
        group = domain.groups[self.ig]
        coors = domain.get_mesh_coors(actual=True)

        if gtype == 'volume':
            qp = self.get_qp('v', integral)

            iels = region.get_cells(self.ig)

            geo_ps = self.interp.get_geom_poly_space('v')
            ps = self.interp.poly_spaces['v']
            bf = self.get_base('v', 0, integral, iels=iels)

            conn = nm.take(group.conn, iels.astype(nm.int32), axis=0)
            mapping = VolumeMapping(coors, conn, poly_space=geo_ps)
            vg = mapping.get_mapping(qp.vals, qp.weights, poly_space=ps,
                                     ori=self.ori)

            out = vg

        elif gtype == 'plate':
            import sfepy.mechanics.membranes as mm
            from sfepy.linalg import dot_sequences

            qp = self.get_qp('v', integral)
            iels = region.get_cells(self.ig)

            ps = self.interp.poly_spaces['v']
            bf = self.get_base('v', 0, integral, iels=iels)

            conn = nm.take(group.conn, nm.int32(iels), axis=0)
            ccoors = coors[conn]

            # Coordinate transformation matrix (transposed!).
            mtx_t = mm.create_transformation_matrix(ccoors)

            # Transform coordinates to the local coordinate system.
            coors_loc = dot_sequences((ccoors - ccoors[:, 0:1, :]), mtx_t)

            # Mapping from transformed elements to reference elements.
            mapping = mm.create_mapping(coors_loc, field.gel, 1)
            vg = mapping.get_mapping(qp.vals, qp.weights, poly_space=ps,
                                     ori=self.ori)
            vg.mtx_t = mtx_t
            out = vg

        elif (gtype == 'surface') or (gtype == 'surface_extra'):
            assert_(field.approx_order > 0)

            if self.ori is not None:
                msg = 'surface integrals do not work yet with the' \
                      ' hierarchical basis!'
                raise ValueError(msg)

            sd = domain.surface_groups[self.ig][region.name]
            esd = self.surface_data[region.name]

            qp = self.get_qp(sd.face_type, integral)

            geo_ps = self.interp.get_geom_poly_space(sd.face_type)
            ps = self.interp.poly_spaces[esd.face_type]
            bf = self.get_base(esd.face_type, 0, integral)

            conn = sd.get_connectivity()

            mapping = SurfaceMapping(coors, conn, poly_space=geo_ps)
            sg = mapping.get_mapping(qp.vals, qp.weights, poly_space=ps,
                                     mode=gtype)
            if gtype == 'surface_extra':
                sg.alloc_extra_data(self.n_ep['v'])

                self.create_bqp(region.name, integral)
                qp = self.qp_coors[(integral.name, esd.bkey)]

                v_geo_ps = self.interp.get_geom_poly_space('v')
                bf_bg = v_geo_ps.eval_base(qp.vals, diff=True)
                ebf_bg = self.get_base(esd.bkey, 1, integral)

                sg.evaluate_bfbgm(bf_bg, ebf_bg, coors, sd.fis, group.conn)

            out =  sg

        elif gtype == 'point':
            out = mapping = None

        else:
            raise ValueError('unknown geometry type: %s' % gtype)

        if out is not None:
            # Store the integral used.
            out.integral = integral
            out.qp = qp
            out.ps = ps
            # Update base.
            out.bf[:] = bf

        if return_mapping:
            out = (out, mapping)

        return out

    def _create_bqp(self, skey, bf_s, weights, integral_name):
        interp = self.interp
        gel = interp.gel
        bkey = 'b%s' % skey[1:]
        bqpkey = (integral_name, bkey)
        coors, faces = gel.coors, gel.get_surface_entities()

        vals = _interp_to_faces(coors, bf_s, faces)
        self.qp_coors[bqpkey] = Struct(name = 'BQP_%s' % bkey,
                                       vals = vals, weights = weights)
        interp.poly_spaces[bkey] = interp.poly_spaces['v']
        return bkey

    def create_bqp(self, region_name, integral):
        sd = self.surface_data[region_name]
        bqpkey = (integral.name, sd.bkey)
        if not bqpkey in self.qp_coors:
            bf_s = self.get_base(sd.face_type, 0, integral,
                                 from_geometry=True)
            qp = self.get_qp(sd.face_type, integral)

            bkey = self._create_bqp(sd.face_type, bf_s, qp.weights,
                                    integral.name)
            assert_(bkey == sd.bkey)

class DiscontinuousApproximation(Approximation):

    def eval_extra_coor(self, coors, mesh_coors):
        """
        Compute coordinates of extra nodes. For discontinuous
        approximations, all nodes are treated as extra.
        """
        gps = self.interp.gel.interp.poly_spaces['v']
        ps = self.interp.poly_spaces['v']

        eval_nodal_coors(coors, mesh_coors, self.region, ps, gps,
                         self.econn, self.ig, only_extra=False)

class SurfaceApproximation(Approximation):

    def __init__(self, name, interp, region, ig):
        Approximation.__init__(self, name, interp, region, ig, is_surface=True)

    def get_qp(self, key, integral):
        """
        Get quadrature points and weights corresponding to the given key
        and integral. The key is 's#', where # is the number of
        face vertices.
        """
        assert_(key[0] == 's')
        qpkey = (integral.name, key)

        if not self.qp_coors.has_key(qpkey):
            interp = self.interp
            geometry = interp.gel.name

            vals, weights = integral.get_qp(geometry)
            self.qp_coors[qpkey] = Struct(vals=vals, weights=weights)

        return self.qp_coors[qpkey]

