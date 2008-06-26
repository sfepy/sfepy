fileName_meshes = ['database/simple.mesh',
                   'database/simple.vtk',
                   'database/t.1.node',
                   'database/maillage.txt']
same = [(0, 1)]

from sfepy.base.testing import TestCommon

##
# c: 05.02.2008
class Test( TestCommon ):
    """Write test names explicitely to impose a given order of evaluation."""
    tests = ['test_readMeshes', 'test_compareSameMeshes']
    
    ##
    # c: 05.02.2008, r: 05.02.2008
    def fromConf( conf, options ):
        return Test( conf = conf, options = options )
    fromConf = staticmethod( fromConf )
    
    ##
    # c: 05.02.2008, r: 05.02.2008
    def test_readMeshes( self ):
        """Try to read all listed meshes."""
        from sfepy.fem.mesh import Mesh

        meshes = {}
        for ii, fileName in enumerate( fileName_meshes ):
            self.report( '%d. mesh: %s' % (ii + 1, fileName) )
            mesh = Mesh.fromFile( fileName )
            self.report( 'read ok' )
            meshes[fileName] = mesh

        self.meshes = meshes

        return True

    ##
    # c: 05.02.2008, r: 05.02.2008
    def test_compareSameMeshes( self ):
        """Compare same meshes in various formats."""
        import numpy as nm

        oks = []
        for i0, i1 in same:
            name0 = fileName_meshes[i0]
            name1 = fileName_meshes[i1]
            self.report( 'comparing meshes from "%s" and "%s"' % (name0, name1) )
            mesh0 = self.meshes[name0]
            mesh1 = self.meshes[name1]

            ok0 = (mesh0.dim == mesh1.dim)
            if not ok0:
                self.report( 'dimension failed!' )
            oks.append( ok0 )
            
            ok0 = mesh0.nNod == mesh1.nNod
            if not ok0:
                self.report( 'number of nodes failed!' )
            oks.append( ok0 )

            ok0 = mesh0.nEl == mesh1.nEl
            if not ok0:
                self.report( 'number of elements failed!' )
            oks.append( ok0 )

            ok0 = mesh0.nEPs == mesh1.nEPs
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

            for ii in range( len( mesh0.matIds ) ):
                ok0 = nm.all( mesh0.matIds[ii] == mesh1.matIds[ii] )
                if not ok0:
                    self.report( 'matIds failed!' )
                oks.append( ok0 )

            for ii in range( len( mesh0.matIds ) ):
                ok0 = nm.all( mesh0.conns[ii] == mesh1.conns[ii] )
                if not ok0:
                    self.report( 'connectivities failed!' )
                oks.append( ok0 )

        return sum( oks ) == len( oks )
