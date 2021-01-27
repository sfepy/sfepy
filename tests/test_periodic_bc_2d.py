# 01.06.2007, c
# last revision: 25.02.2008
from __future__ import absolute_import
from sfepy import data_dir

filename_mesh = data_dir + '/meshes/various_formats/small2d.mesh'

material_1 = {
    'name' : 'coef',
    'values' : {'coef' : 1.0},
}

region_1000 = {
    'name' : 'Omega',
    'select' : 'all',
}
region_1 = {
    'name' : 'Left',
    'select' : 'vertices in (x < -0.499)',
    'kind' : 'facet',
}
region_2 = {
    'name' : 'Right',
    'select' : 'vertices in (x > 0.499)',
    'kind' : 'facet',
}
region_22 = {
    'name' : 'Bottom',
    'select' : 'vertices in (y < -0.499)',
    'kind' : 'facet',
}
region_23 = {
    'name' : 'Top',
    'select' : 'vertices in (y > 0.499)',
    'kind' : 'facet',
}

field_1 = {
    'name' : '2_displacement',
    'dtype' : 'real',
    'shape' : (2,),
    'region' : 'Omega',
    'approx_order' : 2,
}

field_2 = {
    'name' : 'pressure',
    'dtype' : 'real',
    'shape' : (1,),
    'region' : 'Omega',
    'approx_order' : 1,
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

from sfepy.discrete.fem.periodic import match_x_line, match_y_line

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
        from sfepy.discrete import Problem
        problem = Problem.from_conf(conf, init_equations=False)

        test = Test( problem = problem,
                     conf = conf, options = options )
        return test
    from_conf = staticmethod( from_conf )

    ##
    # c: 01.06.2007, r: 18.02.2008
    def test_pbc( self ):
        from sfepy.discrete import Variables, Conditions

        problem  = self.problem
        conf = self.conf

        ebcs = Conditions.from_conf(conf.ebcs, problem.domain.regions)
        epbcs = Conditions.from_conf(conf.epbcs, problem.domain.regions)

        variables = Variables.from_conf(conf.variables, problem.fields)
        variables.equation_mapping(ebcs, epbcs, None, problem.functions)
        state = variables.create_state_vector()
        variables.apply_ebc(state)
        return variables.has_ebc(state)
