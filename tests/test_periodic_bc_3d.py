# 04.06.2007, c
# last revision: 18.02.2008

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

variable_1 = {
    'name' : 'u',
    'kind' : 'unknown field',
    'field' : '3_displacement',
    'dofs' : (0, 1, 2),
    'order' : 0,
}
variable_2 = {
    'name' : 'v',
    'kind' : 'test field',
    'field' : '3_displacement',
    'dofs' : (0, 1, 2),
    'dual' : 'u',
}
variable_3 = {
    'name' : 'p',
    'kind' : 'unknown field',
    'field' : 'pressure',
    'dofs' : (8,),
    'order' : 1,
}
variable_4 = {
    'name' : 'q',
    'kind' : 'test field',
    'field' : 'pressure',
    'dofs' : (8,),
    'dual' : 'p',
}

ebcs = {}
epbc_10 = {
    'name' : 'u_rl',
    'region' : ['Left', 'Right'],
    'dofs' : [(0, 1, 2), (0, 1, 2)],
    'match' : 'matchXPlane',
}
epbc_12 = {
    'name' : 'u_tb',
    'region' : ['Top', 'Bottom'],
    'dofs' : [(0, 1, 2), (0, 1, 2)],
    'match' : 'matchZPlane',
}
epbc_13 = {
    'name' : 'u_nf',
    'region' : ['Near', 'Far'],
    'dofs' : [(0, 1, 2), (0, 1, 2)],
    'match' : 'matchYPlane',
}
epbc_20 = {
    'name' : 'p_rl',
    'region' : ['Left', 'Right'],
    'dofs' : [(8,), (8,)],
    'match' : 'matchXPlane',
}
epbc_22 = {
    'name' : 'p_tb',
    'region' : ['Top', 'Bottom'],
    'dofs' : [(8,), (8,)],
    'match' : 'matchZPlane',
}
epbc_23 = {
    'name' : 'p_nf',
    'region' : ['Near', 'Far'],
    'dofs' : [(8,), (8,)],
    'match' : 'matchYPlane',
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
    # c: 01.06.2007, r: 18.02.2008
    def test_pbc( self ):
        problem  = self.problem
        conf = self.conf
        
        problem.variables.equationMapping( conf.ebcs, conf.epbcs,
                                           problem.domain.regions,
                                           None, conf.funmod )
        state = problem.createStateVector()
        problem.applyEBC( state )
        return problem.variables.hasEBC( state )
