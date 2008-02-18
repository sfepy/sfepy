# 01.06.2007, c
# last revision: 18.02.2008

fileName_mesh = 'database/tests/small2d.mesh'

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
region_22 = {
    'name' : 'Bottom',
    'select' : 'nodes in (y < -0.499)'
}
region_23 = {
    'name' : 'Top',
    'select' : 'nodes in (y > 0.499)'
}

field_1 = {
    'name' : '2_displacement',
    'dim' : (2,1),
    'flags' : (),
    'domain' : 'Omega',
    'bases' : {'Omega' : '2_3_P2'}
}

field_2 = {
    'name' : 'pressure',
    'dim' : (1,1),
    'flags' : (),
    'domain' : 'Omega',
    'bases' : {'Omega' : '2_3_P1'}
}

variable_1 = {
    'name' : 'u',
    'kind' : 'unknown field',
    'field' : '2_displacement',
    'dofs' : (0, 1),
    'order' : 0,
}
variable_2 = {
    'name' : 'v',
    'kind' : 'test field',
    'field' : '2_displacement',
    'dofs' : (0, 1),
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
    'dofs' : [(0, 1), (0, 1)],
    'match' : 'matchYLine',
}
epbc_12 = {
    'name' : 'u_tb',
    'region' : ['Top', 'Bottom'],
    'dofs' : [(0, 1), (0, 1)],
    'match' : 'matchXLine',
}
epbc_20 = {
    'name' : 'p_rl',
    'region' : ['Left', 'Right'],
    'dofs' : [(8,), (8,)],
    'match' : 'matchYLine',
}
epbc_22 = {
    'name' : 'p_tb',
    'region' : ['Top', 'Bottom'],
    'dofs' : [(8,), (8,)],
    'match' : 'matchXLine',
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
