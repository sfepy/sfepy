# 28.05.2009, c
from sfepy import data_dir

filename_mesh = data_dir + '/meshes/3d/elbow.mesh'

fields = {
    'scalar' : ('real', 'scalar', 'Omega', 1),
    'vector' : ('real', 'vector', 'Omega', 1),
}

integrals = {
    'i1' : ('v', 2),
    'i2' : ('s', 2),
}

regions = {
    'Omega' : ('all', {}),
    'Gamma' : ('nodes of surface', {'can_cells' : True}),
}

expressions = {
    'volume_p' : 'd_volume.i1.Omega( p )',
    'volume_u' : 'd_volume.i1.Omega( u )',
    'surface_p' : 'd_volume_surface.i2.Gamma( p )',
    'surface_u' : 'd_volume_surface.i2.Gamma( u )',
}

import numpy as nm
from sfepy.base.testing import TestCommon
from sfepy.base.base import debug, pause

##
# 10.07.2007, c
class Test( TestCommon ):
    tests = ['test_volume']

    def from_conf( conf, options ):
        from sfepy.fem import ProblemDefinition

        problem = ProblemDefinition.from_conf(conf, init_equations=False)
        test = Test( problem = problem,
                     conf = conf, options = options )
        return test
    from_conf = staticmethod( from_conf )

    def test_volume( self ):
        from sfepy.fem import FieldVariable

        ok = True

        field_map = {'u' : 'vector', 'p' : 'scalar'}

        volumes = {}
        avg = 0.0
        for key, term in expressions.items():
            var_name = key[-1]
            field = self.problem.fields[field_map[var_name]]
            var = FieldVariable(var_name, 'parameter', field, 1,
                                primary_var_name='(set-to-None)')

            val = self.problem.evaluate(term, **{var_name : var})

            volumes[key] = val
            avg += val

        avg /= len(volumes)

        for key, val in volumes.items():
            err = nm.abs( avg - val ) / nm.abs( avg )
            _ok = err < 1e-12
            self.report('"'"%s"'" - volume: %e' % (key, val))
            self.report('"'"%s"'" - relative volume difference: %e -> %s'\
                         % (key, err, _ok) )
            ok = ok and _ok

        return ok
