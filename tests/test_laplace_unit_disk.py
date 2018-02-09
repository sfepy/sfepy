# 31.05.2007, c
# last revision: 25.02.2008
from __future__ import absolute_import
from sfepy import data_dir

filename_mesh = data_dir + '/meshes/2d/circle_sym.mesh'

material_1 = {
    'name' : 'coef',
    'values' : {
        'val' : 1.0,
    },
}
material_2 = {
    'name' : 'm',
    'values' : {
        'K' : [[1.0, 0.0], [0.0, 1.0]],
    },
}

field_1 = {
    'name' : 'a_harmonic_field',
    'dtype' : 'real',
    'shape' : 'scalar',
    'region' : 'Omega',
    'approx_order' : 2,
}

variable_1 = {
    'name' : 't',
    'kind' : 'unknown field',
    'field' : 'a_harmonic_field',
    'order' : 0,
}
variable_2 = {
    'name' : 's',
    'kind' : 'test field',
    'field' : 'a_harmonic_field',
    'dual' : 't',
}

region_1000 = {
    'name' : 'Omega',
    'select' : 'all',
}

region_1 = {
    'name' : 'Centre',
    'select' : 'vertices in (x < 1e-8) & (x > -1e-8) & (y < 1e-8) & (y > -1e-8)',
    'kind' : 'vertex'
}

region_2 = {
    'name' : 'Gamma',
    'select' : 'vertices of surface',
    'kind' : 'facet',
}

ebc_1 = {
    'name' : 't_centre',
    'region' : 'Centre',
    'dofs' : {'t.0' : 1.0},
}
ebc_2 = {
    'name' : 't_gamma',
    'region' : 'Gamma',
    'dofs' : {'t.0' : 0.0},
}

integral_1 = {
    'name' : 'i',
    'order' : 2,
}

equations = {
    'Temperature' : """dw_laplace.i.Omega( coef.val, s, t ) = 0"""
}

solution = {
    't' : '- 5.0 * (x - 0.5)',
}

solver_0 = {
    'name' : 'ls',
    'kind' : 'ls.scipy_direct',
}

solver_1 = {
    'name' : 'newton',
    'kind' : 'nls.newton',

    'i_max'      : 1,
    'eps_a'      : 1e-10,
}

from sfepy.base.testing import TestCommon

##
# 31.05.2007, c
class Test( TestCommon ):

    ##
    # 30.05.2007, c
    def from_conf( conf, options ):
        from sfepy.applications import solve_pde

        problem, state = solve_pde(conf, save_results=False)

        test = Test(problem=problem, state=state, conf=conf, options=options)
        return test
    from_conf = staticmethod( from_conf )

    ##
    # 31.05.2007, c
    # 02.10.2007
    def test_boundary_fluxes( self ):
        from sfepy.discrete import Material
        problem = self.problem

        region_names = ['Gamma']

        variables = problem.get_variables()
        get_state = variables.get_state_part_view
        state = self.state.copy(deep=True)

        problem.time_update(ebcs={}, epbcs={})
        ## problem.save_ebc( 'aux.vtk' )

        state.apply_ebc()
        nls = problem.get_nls()
        aux = nls.fun(state())

        field = variables['t'].field

        conf_m = problem.conf.get_item_by_name('materials', 'm')
        m = Material.from_conf(conf_m, problem.functions)

        ok = True
        for ii, region_name in enumerate( region_names ):
            flux_term = 'd_surface_flux.1.%s( m.K, t )' % region_name
            val1 = problem.evaluate(flux_term, t=variables['t'], m=m)

            rvec = get_state( aux, 't', True )
            reg = problem.domain.regions[region_name]
            nods = field.get_dofs_in_region(reg, merge=True)
            val2 = rvec[nods].sum() # Assume 1 dof per node.

            eps = 1e-2
            ok = ok and ((abs( val1 - val2 ) < eps))
            self.report( '%d. %s: |%e - %e| = %e < %.2e'\
                         % (ii, region_name, val1, val2, abs( val1 - val2 ),
                            eps) )

        return ok
