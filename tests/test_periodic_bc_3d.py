# 04.06.2007, c
# last revision: 10.12.2007

fileName_mesh = 'database/tests/small3d.mesh'

material_1 = {
    'name' : 'coef',
    'mode' : 'here',
    'region' : 'Omega',
    'coef' : 1.0,
}

region_1000 = {
    'name' : 'Omega',
    'select' : 'all',
}
region_1 = {
    'name' : 'Left',
    'select' : 'nodes in (x < -0.499)',
}
region_2 = {
    'name' : 'Right',
    'select' : 'nodes in (x > 0.499)',
}
region_3 = {
    'name' : 'Near',
    'select' : 'nodes in (y < -0.499)',
}
region_4 = {
    'name' : 'Far',
    'select' : 'nodes in (y > 0.499)',
}
region_5 = {
    'name' : 'Bottom',
    'select' : 'nodes in (z < -0.499)'
}
region_6 = {
    'name' : 'Top',
    'select' : 'nodes in (z > 0.499)'
}

field_1 = {
    'name' : '3_displacement',
    'dim' : (3,1),
    'flags' : (),
    'domain' : 'Omega',
    'bases' : {'Omega' : '3_4_P2'}
}

field_2 = {
    'name' : 'pressure',
    'dim' : (1,1),
    'flags' : (),
    'domain' : 'Omega',
    'bases' : {'Omega' : '3_4_P1'}
}

variables = {
    'u'   : ('field', 'unknown', '3_displacement', (0, 1, 2), 0),
    'v'   : ('field', 'test',    '3_displacement', (0, 1, 2), 'u'),
    'p'   : ('field', 'unknown', 'pressure',       (8,), 1),
    'q'   : ('field', 'test',    'pressure',       (8,), 'p'),
}

ebc = {}
epbc = {
    'Left' : (('u', (0, 1, 2), (0, 1, 2), 'Right', 'matchXPlane'),
              ('p', (8,), (8,), 'Right', 'matchXPlane')),
    'Near' : (('u', (0, 1, 2), (0, 1, 2), 'Far', 'matchYPlane'),
              ('p', (8,), (8,), 'Far', 'matchYPlane')),
    'Top'  : (('u', (0, 1, 2), (0, 1, 2), 'Bottom', 'matchZPlane'),
              ('p', (8,), (8,), 'Bottom', 'matchZPlane')),
}

fe = {
    'chunkSize' : 1000
}

from sfe.fem.periodic import *
from sfe.base.testing import TestCommon

##
# 01.06.2007, c
class Test( TestCommon ):

    ##
    # 01.06.2007, c
    def fromConf( conf, options ):
        from sfe.fem.problemDef import ProblemDefinition
        problem = ProblemDefinition.fromConf( conf, initEquations = False )

        test = Test( problem = problem,
                     conf = conf, options = options )
        return test
    fromConf = staticmethod( fromConf )

    ##
    # 01.06.2007, c
    def test_pbc( self ):
        problem  = self.problem
        conf = self.conf
        
        problem.variables.equationMapping( conf.ebc, conf.epbc,
                                           problem.domain.regions,
                                           None, conf.funmod )
        state = problem.createStateVector()
        problem.applyEBC( state )
        return problem.variables.hasEBC( state )
