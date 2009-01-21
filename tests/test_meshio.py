filename_meshes = ['database/simple.mesh',
                   'database/simple.vtk',
                   'database/t.1.node',
                   'database/maillage.txt',
                   'database/tests/hex4.mesh3d',
                   'database/tests/tetra8.mesh3d']
same = [(0, 1)]

from sfepy.base.testing import TestCommon, assert_

##
# c: 05.02.2008
class Test( TestCommon ):
    """Write test names explicitely to impose a given order of evaluation."""
    tests = ['test_read_meshes', 'test_compare_same_meshes', 'test_read_dimension']
    
    ##
    # c: 05.02.2008, r: 05.02.2008
    def from_conf( conf, options ):
        return Test( conf = conf, options = options )
    from_conf = staticmethod( from_conf )
    
    ##
    # c: 05.02.2008, r: 05.02.2008
    def test_read_meshes( self ):
        """Try to read all listed meshes."""
        from sfepy.fem import Mesh

        meshes = {}
        for ii, filename in enumerate( filename_meshes ):
            self.report( '%d. mesh: %s' % (ii + 1, filename) )
            mesh = Mesh.from_file( filename )

            assert_(mesh.dim == (mesh.nod0.shape[1]-1))
            assert_(mesh.n_nod == (mesh.nod0.shape[0]))
            assert_(mesh.n_el == sum(mesh.n_els))
            for ig, conn in enumerate( mesh.conns ):
                assert_(conn.shape[0] == len(mesh.mat_ids[ig]))
                assert_(conn.shape[0] == mesh.n_els[ig])
                assert_(conn.shape[1] == mesh.n_e_ps[ig])
                
            self.report( 'read ok' )
            meshes[filename] = mesh

        self.meshes = meshes

        return True

    ##
    # c: 05.02.2008, r: 05.02.2008
    def test_compare_same_meshes( self ):
        """Compare same meshes in various formats."""
        import numpy as nm

        oks = []
        for i0, i1 in same:
            name0 = filename_meshes[i0]
            name1 = filename_meshes[i1]
            self.report( 'comparing meshes from "%s" and "%s"' % (name0, name1) )
            mesh0 = self.meshes[name0]
            mesh1 = self.meshes[name1]

            ok0 = (mesh0.dim == mesh1.dim)
            if not ok0:
                self.report( 'dimension failed!' )
            oks.append( ok0 )
            
            ok0 = mesh0.n_nod == mesh1.n_nod
            if not ok0:
                self.report( 'number of nodes failed!' )
            oks.append( ok0 )

            ok0 = mesh0.n_el == mesh1.n_el
            if not ok0:
                self.report( 'number of elements failed!' )
            oks.append( ok0 )

            ok0 = mesh0.n_e_ps == mesh1.n_e_ps
            if not ok0:
                self.report( 'number of element points failed!' )
            oks.append( ok0 )

            ok0 = mesh0.descs == mesh1.descs
            if not ok0:
                self.report( 'element types failed!' )
            oks.append( ok0 )

            ok0 = nm.allclose( mesh0.nod0, mesh1.nod0 )
            if not ok0:
                self.report( 'nodes failed!' )
            oks.append( ok0 )

            for ii in range( len( mesh0.mat_ids ) ):
                ok0 = nm.all( mesh0.mat_ids[ii] == mesh1.mat_ids[ii] )
                if not ok0:
                    self.report( 'mat_ids failed!' )
                oks.append( ok0 )

            for ii in range( len( mesh0.mat_ids ) ):
                ok0 = nm.all( mesh0.conns[ii] == mesh1.conns[ii] )
                if not ok0:
                    self.report( 'connectivities failed!' )
                oks.append( ok0 )

        return sum( oks ) == len( oks )

    ##
    # c: 03.07.2008, r: 03.07.2008
    def test_read_dimension( self ):
        from sfepy.fem import MeshIO
        meshes = {'database/tests/small2d.mesh' : 2,
                  'database/tests/small2d.vtk' : 2,
                  'database/tests/small3d.mesh' : 3,
                  'database/simple.mesh' : 3,
                  'database/simple.vtk' : 3}

        ok = True
        for filename, adim in meshes.iteritems():
            self.report( 'mesh: %s, dimension %d' % (filename, adim) )
            io = MeshIO.any_from_filename( filename )
            dim = io.read_dimension()
            if dim != adim:
                self.report( 'read dimension %d -> failed' % dim )
                ok = False
            else:
                self.report( 'read dimension %d -> ok' % dim )

        return ok
