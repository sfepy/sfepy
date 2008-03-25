"""
This module is not a test file. It contains classes grouping some common
functionality, that is used in several test files.
"""

from sfe.base.testing import TestCommon
import os.path as op

##
# 05.06.2007, c
class TestInput( TestCommon ):
    """Test that an input file works. See test_input_*.py files."""

    ##
    # c: 05.06.2007, r: 19.02.2008
    def fromConf( conf, options, cls = None ):
        from sfe.base.conf import ProblemConf, getStandardKeywords

        required, other = getStandardKeywords()
        testConf = ProblemConf.fromFile( conf.inputName, required, other )

        if cls is None:
            cls = TestInput
        test = cls( testConf = testConf, conf = conf, options = options )

        return test
    fromConf = staticmethod( fromConf )

    ##
    # c: 05.06.2007, r: 08.02.2008
    def test_input( self ):
        from sfe.solvers.generic import solveStationary
        from sfe.base.base import IndexedStruct

        self.report( 'solving %s...' % self.conf.inputName )
        status = IndexedStruct()
        dpb, vecDP, data = solveStationary( self.testConf, nlsStatus = status )
        ok = status.condition == 0
        out = dpb.stateToOutput( vecDP )

        name = op.join( self.options.outDir,
                        op.split( self.conf.outputName )[1] )
        dpb.saveState( name, vecDP )
        self.report( '%s solved' % self.conf.inputName )

        return ok

##
# 03.10.2007, c
class TestLCBC( TestCommon ):
    """Test linear combination BC. See test_lcbc_*.py files."""

    ##
    # 03.10.2007, c
    def fromConf( conf, options ):
        return TestLCBC( conf = conf, options = options )
    fromConf = staticmethod( fromConf )

    ##
    # c: 03.10.2007, r: 08.02.2008
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

        status = IndexedStruct()
        problem, vec, data = solveStationary( self.conf,
                                              nlsStatus = status )
        ok = status.condition == 0
        self.report( 'converged: %s' % ok )
        out = problem.stateToOutput( vec )

        strain = evalTermOP( vec, 'de_cauchy_strain.i1.Y( u )', problem )
        out['strain'] = Struct( name = 'output_data',
                                mode = 'cell', data = strain,
                                dofTypes = None )

        name = op.join( self.options.outDir,
                        op.split( self.conf.outputName )[1] )
        problem.domain.mesh.write( name, io = 'auto', out = out )

        ##
        # Check if rigid body displacements are really rigid should go here.

        return ok
