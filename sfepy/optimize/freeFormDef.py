import tables as pt
from sfepy.base.base import *
from sfepy.linalg.utils import cycle
from sfepy.fem.mesh import Mesh

##
# 11.01.2006, c
# 20.01.2006
# 23.01.2006
# 24.01.2006
# 12.04.2006
def read_spline_box_hdf5( filename ):
    if not pt.isHDF5File( filename ):
        raise ValueError, 'not a HDF5 file! (%s)' % filename

    fd = pt.openFile( filename, mode = 'r' )
    boxes = fd.listNodes( '/box' )
    n_box = len( boxes )
    dim = len( fd.listNodes( boxes[0].ax ) )

    sp_boxes = SplineBoxes( dim = dim, n_box = n_box, n_vertex = 0,
                           spbs = OneTypeList( SplineBox ) )
    for box in boxes:
        spb = SplineBox()
        sp_boxes.spbs.append( spb )

        spb.ib = int( box._v_name )
        spb.cpi = nm.asarray( box.cpi.read() ) - 1
        spb.gpi = nm.asarray( box.gpi.read() ) - 1
        spb.cxyz = nm.asarray( box.cxyz.read() ).transpose()
        spb.cxyz0 = spb.cxyz.copy()
        spb.ax = []
        for axi in fd.listNodes( box.ax ):
            spb.ax.append( nm.asarray( axi.bsc.read() ) )

        sp_boxes.n_vertex = max( sp_boxes.n_vertex, nm.amax( spb.gpi ) + 1 )
        print nm.amin( spb.gpi ), nm.amax( spb.gpi )

        ##
        # Fix cpi by rebuilding :).
        off = 0
        n0, n1, n2 = spb.cpi.shape
        aux = nm.arange( n0 * n1 ).reshape( n1, n0 ).transpose()
        for ii in xrange( n2 ):
            spb.cpi[:,:,ii] = aux + off
            off += n0 * n1

    fd.close()

    for perm in cycle( [n_box] * 2 ):
        if perm[0] == perm[1]: continue
        gpi1 = sp_boxes.spbs[perm[0]].gpi
        gpi2 = sp_boxes.spbs[perm[1]].gpi
        assert_( len( nm.intersect1d( gpi1, gpi2 ) ) == 0 )
    
    return sp_boxes

##
# 20.01.2006, c
# 17.02.2006
# 12.04.2006
# 13.04.2006
def read_dsg_vars_hdf5( filename ):
    if not pt.isHDF5File( filename ):
        raise ValueError, 'not a HDF5 file! (%s)' % filename

    fd = pt.openFile( filename, mode = 'r' )
    aux1 = fd.getNode( '/dsgvar/inx' ).read()
    aux2 = fd.getNode( '/dsgvar/val' ).read()
    aux3 = fd.getNode( '/dsgvar/nsbdsg' ).read()
    dsg_vars = DesignVariables( indx = nm.asarray( aux1, dtype = nm.int32 ),
                               cxyz =  nm.asarray( aux2 ),
                               null_space_b = nm.asarray( aux3 ) )
    dsg_vars.indx = dsg_vars.indx.transpose()
    dsg_vars.indx[:,1] -= 1

    # No. of design variables.
    dsg_vars.n_dsg = dsg_vars.null_space_b.shape[1]
    # No. of control points.
    dsg_vars.n_cp = dsg_vars.indx.shape[0]
    # Design vector and initial design vector. 
    dsg_vars.val0 = nm.zeros( (dsg_vars.n_dsg,), dtype = nm.float64 )
    dsg_vars.val = dsg_vars.val0.copy()
    
    fd.close()

    return dsg_vars

##
# 12.04.2006, c
def interp_box_coordinates( spb, cxyz = None ):
    if cxyz is None:
        cxyz = spb.cxyz
        
    dim = len( spb.ax )
    pp = nm.zeros( (spb.ax[0].shape[0], dim), dtype = nm.float64 )
    for ii in xrange( spb.cpi.shape[0] ):
        for ij in xrange( spb.cpi.shape[1] ):
            aux = spb.ax[0][:,ii] * spb.ax[1][:,ij]
            for ik in xrange( spb.cpi.shape[2] ):
                aux2 = aux * spb.ax[2][:,ik]
                ip = spb.cpi[ii,ij,ik]
                pp += aux2[:,nm.newaxis] * cxyz[ip,:]
##                         print ii, ij, ik, ip, cxyz[:,ip]
##                         pause()
    return pp

class DesignVariables( Struct ):
    ##
    # 20.01.2006, c
    def renumber_by_boxes( self, sp_boxes ):
        bmap = nm.array( [[ii, spb.ib + 1]
                          for ii, spb in enumerate( sp_boxes.spbs )],
                         dtype = nm.int32 )
        bpos = nm.zeros( (nm.amax( bmap[:,1] ) + 1,), dtype = nm.int32 )
        bpos[bmap[:,1]] = bmap[:,0]
        print bmap, bpos
        self.indx[:,0] = bpos[self.indx[:,0]]

    ##
    # 24.04.2006, c
    def normalize_null_space_base( self, magnitude = 1.0 ):
        for ic in xrange( self.n_dsg ):
            col_norm = nla.norm( self.null_space_b[:,ic] )
            self.null_space_b[:,ic] /= magnitude * col_norm

class SplineBox( Struct ):
    pass

class SplineBoxes( Struct ):
    ##
    # 11.01.2006, c
    # 20.01.2006
    # 12.04.2006
    def interp_mesh_velocity( self, shape, dsg_vars, idsg ):
        n_nod, dim = shape
        cv = dsg_vars.null_space_b[:,idsg]
        cv = nm.reshape( cv, (dsg_vars.n_cp, dim) )

        vel = nm.zeros( shape = shape, dtype = nm.float64 )

        for ib, spb in enumerate( self.spbs ):
            ii = nm.where( ib == dsg_vars.indx[:,0] )[0]
#            print ib, ii
            if ii.shape[0] > 0:
                delta = nm.zeros_like( spb.cxyz )
                delta[dsg_vars.indx[ii,1],:] = cv[ii,:]

                vv = interp_box_coordinates( spb, delta )
##                 print delta
##                 print vv
##                 pause()
                vel[spb.gpi,:] = vv
        return vel

    ##
    # 23.01.2006, c
    # 24.01.2006
    # 12.04.2006
    def interp_coordinates( self ):
        coors = nm.zeros( (self.n_vertex, self.dim), nm.float64 )
        # Body using the (design-deformed) Greville abscisae spb.cxyz.
        for spb in self.spbs:
            coors[spb.gpi] = interp_box_coordinates( spb )

        return coors

    ##
    # 13.04.2006, c
    # 19.04.2006
    def set_control_points( self, dsg_vars = None ):
        if dsg_vars is None:
            # Restore Greville abscisae spb.cxyz0.
            for spb in self.spbs:
                spb.cxyz = spb.cxyz0.copy()
        else:
            # Set design-deformed Greville abscisae spb.cxyz.
            delta = nm.dot( dsg_vars.null_space_b, dsg_vars.val )
            delta = nm.reshape( delta, (dsg_vars.n_cp, self.dim) )
            for ib, spb in enumerate( self.spbs ):
                spb.cxyz = spb.cxyz0.copy()

                ii = nm.where( ib == dsg_vars.indx[:,0] )[0]
                if ii.shape[0] > 0:
                    ir = dsg_vars.indx[ii,1]
                    spb.cxyz[ir,:] = spb.cxyz0[ir,:] + delta[ii,:]

    ##
    # 27.03.2007, c
    def create_mesh_from_control_points( self ):
        offset = 0
        dim = self.spbs[0].cxyz.shape[1]
        coors = nm.empty((0, dim), dtype=nm.float64)
        conns = []
        mat_ids = []
        descs = []
        for ib, spb in enumerate( self.spbs ):
            n_nod = spb.cxyz.shape[0]
            coors = nm.concatenate( (coors, spb.cxyz), 0 )
            descs.append( '3_2' )

            conn = []
            for ij in xrange( spb.cpi.shape[1] ):
                for ik in xrange( spb.cpi.shape[2] ):
                    inx = spb.cpi[:,ij,ik]
                    row = [[p1, p2] for p1, p2 in zip( inx[:-1], inx[1:] )]
                    conn.extend( row )
            for ij in xrange( spb.cpi.shape[0] ):
                for ik in xrange( spb.cpi.shape[2] ):
                    inx = spb.cpi[ij,:,ik]
                    row = [[p1, p2] for p1, p2 in zip( inx[:-1], inx[1:] )]
                    conn.extend( row )
            for ij in xrange( spb.cpi.shape[0] ):
                for ik in xrange( spb.cpi.shape[1] ):
                    inx = spb.cpi[ij,ik,:]
                    row = [[p1, p2] for p1, p2 in zip( inx[:-1], inx[1:] )]
                    conn.extend( row )

            aux = nm.empty(len(conn), dtype=nm.int32)
            aux.fill(ib)
            mat_ids.append(aux)

            conns.append( offset + nm.array( conn, dtype = nm.int32 ) )
            offset += n_nod

        mesh = Mesh.from_data('control_points', coors, None, conns,
                              mat_ids, descs)
        return mesh

if __name__ == '__main__':
    filename = 'klik.h5'
    spboxes = read_spline_box_hdf5( filename )
    dsg_vars = read_dsg_vars_hdf5( filename )
