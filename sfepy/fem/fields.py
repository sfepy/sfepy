from sfepy.base.base import *
import sfepy.base.la as la
import fea
from mesh import Mesh, make_point_cells
import sfepy.terms as terms
import extmods.geometry as gm

##
# 14.07.2006, c
class Fields( Container ):

    ##
    # 14.07.2006, c
    def from_conf(conf, regions):

        objs = OneTypeList( Field )
        for key, val in conf.iteritems():
            objs.append(Field.from_conf(val, regions))

        obj = Fields( objs )

        return obj
    from_conf = staticmethod( from_conf )

    ##
    # 19.07.2006, c
    def setup_coors( self ):
        for field in self:
            field.setup_coors()

##
# 14.07.2006, c
class Field( Struct ):
    all_interps = {}

    def from_conf(conf, regions):
        """To refactor... very hackish now."""
        approx_order = conf.bases.values()[0][5:]
        obj = Field(name = conf.name,
                    dtype = getattr( conf, 'dtype', nm.float64 ),
                    shape = conf.dim[:1],
                    region = regions[conf.domain],
                    approx_order = approx_order)
        return obj
    from_conf = staticmethod( from_conf )

    def __init__(self, name, dtype, shape, region,
                 space='H1', poly_space_base='lagrange', approx_order=1):
        """Create a Field.

        Parameters
        ----------
        name : str
            Object name.
        dtype : numpy.dtype
            Field data type: float64 or complex128.
        shape : int/str
            Field shape: 1 or 'scalar', space dimension (2 or 3) or 'vector'.
            The field shape determines the shape of the FE base functions and
            can be different from a FieldVariable instance shape. 
        region : Region
            The region where the field is defined.
        space : str
            The function space name.
        poly_space_base : str
            The name of polynomial space base.
        approx_order : int/str
            FE approxiamtion order, e.g. 0, 1, 2, '1B' (1 with bubble).
        """
        if isinstance(shape, str):
            try:
                shape = {'scalar' : 1,
                         'vector' : region.domain.shape.dim}[shape]
            except KeyError:
                raise ValueError('unsupported field shape! (%s)', shape)

        Struct.__init__(self,
                        name = name,
                        dtype = dtype,
                        shape = shape,
                        region = region,
                        space = space,
                        poly_space_base = poly_space_base)
        self.domain = self.region.domain

        self.set_approx_order(approx_order)

        # To refactor below...
        self.setup_bases()
        self.create_interpolants()
        self.setup_approximations()
##         print self.aps
##         pause()
        self.setup_global_base()
        self.setup_coors()

    def set_approx_order(self, approx_order):
        """Set a uniform approximation order."""
        
        ao_msg = 'unsupported approximation order! (%s)'
        force_bubble = False

        try:
            ao = int(approx_order)
        except ValueError:
            if approx_order[-1] == 'B':
                ao = int(approx_order[:-1])
                force_bubble = True
            else:
                raise ValueError(ao_msg % approx_order)

        if ao < 0:
            raise ValueError(ao_msg % approx_order)

        self.approx_order = '%s' % approx_order
        self._ao = ao
        self.force_bubble = force_bubble

    def setup_bases(self):
        """Setup FE bases according to self.approx_order and region cell
        types. Assumes one cell type for the whole region!"""
        dim = self.domain.shape.dim
        group = self.domain.groups[self.region.igs[0]]
        n_ep = group.shape.n_ep

        if n_ep == (dim + 1): # simplex
            kind = 'P'
        else: # tensor product
            kind = 'Q'

        self.bases = {self.region.name :
                      '%d_%d_%s%s' % (dim, n_ep, kind, self.approx_order)}

    def create_interpolants(self):
        """
        One interpolant for each base, shared if same.
        """
        self.interps = {}
        for region_name, base_name in self.bases.iteritems():

            if base_name in self.all_interps:
                interp = self.all_interps[base_name]

            else:
                gel = self.domain.geom_els[base_name[:3]]
                interp = fea.Interpolant(base_name, gel)
                self.all_interps[base_name] = interp

            self.interps[base_name] = interp
            
    ##
    # 18.07.2006, c
    # 05.09.2006
    def setup_approximations(self):
        self.aps = fea.Approximations(self.bases, self.interps, self.domain)

    ##
    #
    def igs( self ):
        return self.aps.igs

    ##
    # 19.07.2006, c
    def setup_global_base( self ):

        self.aps.describe_nodes()
        self.aps.setup_nodes()

        aux = self.aps.setup_global_base( self.domain )
        self.n_nod, self.remap, self.cnt_vn, self.cnt_en = aux

##         print self.n_nod, self.cnt_vn, self.cnt_en
#        pause()

    ##
    # 19.07.2006, c
    def setup_coors( self ):
        """Coordinates of field nodes."""
        self.aps.setup_coors( self.domain.mesh, self.cnt_vn )

    ##
    # c: 02.01.2008, r: 02.01.2008
    def get_extra_nodes_as_simplices( self, iextra = None ):
        dim = self.domain.mesh.dim
        if iextra is None:
            noft = self.aps.node_offset_table
            iextra = nm.arange( noft[1,0], noft[-1,-1], dtype = nm.int32 )
        extra = make_point_cells( iextra, dim )
        return {2 : '2_3', 3 : '3_4'}[dim], -nm.ones_like( iextra ), extra

    def create_mesh(self, extra_nodes=True):
        """
        Create a mesh from the field region, optionally including the field
        extra nodes.
        """
        mesh = self.domain.mesh

        if self.approx_order != '0':
            conns, mat_ids, descs = [], [], []
            for region_name, ig, ap in self.aps.iter_aps():
                region = ap.region
                group = region.domain.groups[ig]
                if extra_nodes:
                    conn = ap.econn
                else:
                    offset = group.shape.n_ep
                    conn = ap.econn[:,:offset]
                conns.append(conn)
                mat_ids.append(mesh.mat_ids[ig])
                descs.append(mesh.descs[ig])

            mesh = Mesh.from_data(self.name,
                                  self.aps.coors, None, conns,
                                  mat_ids, descs)

        return mesh

    ##
    # c: 19.07.2006, r: 27.02.2008
    def write_mesh( self, name_template, field_name = None ):
        """Extra nodes are written as zero-size simplices (= points)."""
        if field_name is None:
            field_name = self.name

        tmp = self.create_mesh(extra_nodes=False)

        aux = self.get_extra_nodes_as_simplices()
        tmp.descs.append( aux[0] )
        tmp.mat_ids.append( aux[1] )
        tmp.conns.append( aux[2] )

##         print tmp
##         pause()
        tmp.write( io = 'auto' )

    ##
    # c: 20.07.2006, r: 15.01.2008
    def get_node_descs( self, region ):
        nds = {}
        for region_name, ig, ap in self.aps.iter_aps():
            if ig in region.igs:
                nds[ig] = self.aps.node_descs[region_name]

        return nds

    ##
    # Modify me for bubble-only approximations to not generate vertex nodes.
    # 12.10.2005, c
    # 26.10.2005
    # 26.05.2006
    # 05.06.2006
    # 25.07.2006
    # 04.09.2006
    def interp_c_vals_to_n_vals( self, vec ):
        """len( vec ) == domain.n_el"""
        n_els = [sub.n_el for sub in self.domain.subs]
        oel = nm.cumsum( [0] + n_els )
        if sum( n_els ) != vec.shape[0]:
            print 'incomatible shape! (%d == %d)' % (sum( n_els ), vec.shape[0])
            raise ValueError

        ##
        # Mesh vertex values. 
        n_vertex = self.domain.n_nod
        nod_vol = nm.zeros( (n_vertex,), nm.float64 )
        dim = vec.shape[1]
        nod_vol_val = nm.zeros( (n_vertex, dim ), nm.float64 )
        for ii, ap in enumerate( self.aps ):
            sub = ap.sub
            ig = sub.iseq
            vg = self.vgs[ii]
            volume = nm.squeeze( vg.variable( 2 ) );

            for ii in range( sub.conn.shape[1] ):
                cc = sub.conn[:,ii]
                nod_vol[cc] += volume
                val = volume[:,nm.newaxis] * vec[oel[ig]:oel[ig+1],:]
                ind2, ind1 = nm.meshgrid( nm.arange( dim ), cc )
                nod_vol_val[ind1,ind2] += val
        
        nod_vol_val = nod_vol_val / nod_vol[:,nm.newaxis]

        ##
        # Field nodes values.
        enod_vol_val = self.interp_v_vals_to_n_vals( nod_vol_val )

        return enod_vol_val
    
    ##
    # 05.06.2006, c
    # 25.07.2006
    # 31.08.2006
    def interp_v_vals_to_n_vals( self, vec ):
        dim = vec.shape[1]
        enod_vol_val = nm.zeros( (self.n_nod, dim), nm.float64 )
        for region_name, ig, ap in self.aps.iter_aps():
            group = self.domain.groups[ig]
            offset = group.shape.n_ep
            conn = ap.econn[:,:offset]

            noff = ap.node_offsets.ravel()
            if noff[1] == noff[-1]:
                # Vertex values only...
                ii = nm.unique1d(conn) # Probably wrong?!
                enod_vol_val[ii] = vec[ii]
                continue

            econn = ap.econn

            ginterp = ap.interp.gel.interp
            coors = ap.interp.poly_spaces['v'].node_coors

            bf = ginterp.poly_spaces['v'].eval_base(coors)
            bf = bf[:,0,:].copy()
            
            fea.mu.interp_vertex_data(enod_vol_val, econn, vec, group.conn,
                                      bf, 0)

        return enod_vol_val

    ##
    # 08.08.2006, c
    # 13.02.2007
    def get_coor( self, nods = None, igs = None ):
        """Will igs be ever needed?"""
        if nods is None:
            return self.aps.coors
        else:
            return self.aps.coors[nods]
