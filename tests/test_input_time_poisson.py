# c: 06.02.2008, r: 06.02.2008
inputName = 'input/time_poisson.py'
outputNameTrunk = 'test_time_poisson'

from testsBasic import TestInput

##
# c: 06.02.2008
class Test( TestInput ):

    ##
    # c: 06.02.2008, r: 06.02.2008
    def fromConf( conf, options ):
        return TestInput.fromConf( conf, options, cls = Test )
    fromConf = staticmethod( fromConf )

    ##
    # c: 06.02.2008, r: 06.02.2008
    def test_input( self ):
        """Does not verify anything!!!"""
        import os.path as op
        from sfe.solvers.generic import solveDirect
        from sfe.base.base import Struct

        self.report( 'solving %s...' % self.conf.inputName )
        name = op.join( self.options.outDir, self.conf.outputNameTrunk )
        options = Struct( outputFileNameTrunk = name,
                          saveEBC = False, saveRegions = False,
                          saveFieldMeshes = False, saveRegionFieldMeshes = False,
                          solveNot = False )
        dpb, vecDP, data = solveDirect( self.testConf, options )

        self.report( '%s solved' % self.conf.inputName )

        return True
