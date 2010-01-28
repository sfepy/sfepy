# 01.06.2007, c
# last revision: 25.02.2008
from sfepy import top_dir

filename_mesh = top_dir + '/meshes/various_formats/small2d.mesh'

material_1 = {
    'name' : 'coef',
    'region' : 'Omega',
    'values' : {'coef' : 1.0},
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
    'order' : 0,
}
variable_2 = {
    'name' : 'v',
    'kind' : 'test field',
    'field' : '2_displacement',
    'dual' : 'u',
}
variable_3 = {
    'name' : 'p',
    'kind' : 'unknown field',
    'field' : 'pressure',
    'order' : 1,
}
variable_4 = {
    'name' : 'q',
    'kind' : 'test field',
    'field' : 'pressure',
    'dual' : 'p',
}

ebcs = {}
epbc_10 = {
    'name' : 'rl',
    'region' : ['Left', 'Right'],
    'dofs' : {'u.all' : 'u.all', 'p.0' : 'p.0'},
    'match' : 'match_y_line',
}
epbc_12 = {
    'name' : 'tb',
    'region' : ['Top', 'Bottom'],
    'dofs' : {'u.all' : 'u.all', 'p.0' : 'p.0'},
    'match' : 'match_x_line',
}

fe = {
    'chunk_size' : 1000
}

from sfepy.fem.periodic import *

functions = {
    'match_x_line' : (match_x_line,),
    'match_y_line' : (match_y_line,),
}

from sfepy.base.testing import TestCommon

##
# 01.06.2007, c
class Test( TestCommon ):

    ##
    # 01.06.2007, c
    def from_conf( conf, options ):
        from sfepy.fem import ProblemDefinition
        problem = ProblemDefinition.from_conf( conf, init_equations = False )

        test = Test( problem = problem,
                     conf = conf, options = options )
        return test
    from_conf = staticmethod( from_conf )

    ##
    # c: 01.06.2007, r: 18.02.2008
    def test_pbc( self ):
        problem  = self.problem
        conf = self.conf
        
        problem.variables.equation_mapping(conf.ebcs, conf.epbcs,
                                           problem.domain.regions,
                                           None, problem.functions)
        state = problem.create_state_vector()
        problem.apply_ebc( state )
        return problem.variables.has_ebc( state )
