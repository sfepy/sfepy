from sfe.base.testing import TestCommon
import os.path as op

##
# 05.06.2007, c
class TestInput( TestCommon ):

    ##
    # 05.06.2007, c
    def fromConf( conf, options ):
        from sfe.base.conf import ProblemConf

        required = ['fileName_mesh', 'field_[0-9]+', 'ebc|nbc', 'fe',
                    'equations', 'region_[0-9]+', 'variables',
                    'material_[0-9]+', 'solver_[0-9]+']
        other = ['functions', 'modules', 'epbc', 'lcbc']
        testConf = ProblemConf.fromFile( conf.inputName, required, other )
        test = TestInput( testConf = testConf, conf = conf, options = options )
        return test
    fromConf = staticmethod( fromConf )

    ##
    # 05.06.2007, c
    # 03.07.2007
    # 19.07.2007
    def test_input( self ):
        from sfe.solvers.generic import solveStationary
        from sfe.base.base import IndexedStruct
        import sfe.base.ioutils as io

        self.report( 'solving %s...' % self.conf.inputName )
        status = IndexedStruct()
        dpb, vecDP, data = solveStationary( self.testConf, nlsStatus = status )
        ok = status.condition == 0
        out = dpb.stateToOutput( vecDP )

        name = op.join( self.options.outDir,
                        op.split( self.conf.outputName )[1] )
        fd = open( name, 'w' )
        io.writeVTK( fd, dpb.domain.mesh, out )
        fd.close()
        self.report( '%s solved' % self.conf.inputName )

        return ok

##
# 03.10.2007, c
class TestLCBC( TestCommon ):

    ##
    # 03.10.2007, c
    def fromConf( conf, options ):
        return TestLCBC( conf = conf, options = options )
    fromConf = staticmethod( fromConf )

    ##
    # 03.10.2007, c
    # 05.10.2007
    def test_linearRigidBodyBC( self ):
        import scipy
        if scipy.version.version == "0.6.0":
            # This test uses a functionality implemented in scipy svn, which is
            # missing in scipy 0.6.0
            return True
        from sfe.base.base import Struct
        from sfe.solvers.generic import solveStationary
        from sfe.base.base import IndexedStruct
        from sfe.fem.evaluate import evalTermOP
        import sfe.base.ioutils as io

        status = IndexedStruct()
        problem, vec, data = solveStationary( self.conf,
                                              nlsStatus = status )
        ok = status.condition == 0
        self.report( 'converged: %s' % ok )
        out = problem.stateToOutput( vec )

        strain = evalTermOP( vec, 'de_sdcc_strain.i1.Y( u )', problem )
        out['strain'] = Struct( name = 'output_data',
                                mode = 'cell', data = strain,
                                dofTypes = None )

        name = op.join( self.options.outDir,
                        op.split( self.conf.outputName )[1] )
        fd = open( name, 'w' )
        io.writeVTK( fd, problem.domain.mesh, out )
        fd.close()

        ##
        # Check if rigid body displacements are really rigid should go here.

        return ok
