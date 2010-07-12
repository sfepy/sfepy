# 28.05.2009, c
from sfepy import data_dir

filename_mesh = data_dir + '/meshes/3d/elbow.mesh'

fields = {
    'scalar' : ((1,1), 'real', 'Omega', {'Omega' : '3_4_P1'}),
    'vector' : ((3,1), 'real', 'Omega', {'Omega' : '3_4_P1'}),
}

variables = {
    'p' : ('parameter field', 'scalar', None),
    'u' : ('parameter field', 'vector', None),
}

integrals = {
    'i1' : ('v', 'gauss_o2_d3'),
    'i2' : ('s', 'gauss_o2_d2'),
}

regions = {
    'Omega' : ('all', {}),
    'Gamma' : ('nodes of surface', {'can_cells' : True}),
}

expressions = {
    'volume_s': 'd_volume.i1.Omega( p )',
    'volume_u': 'd_volume.i1.Omega( u )',
    'surface' : 'd_volume_surface.i2.Gamma( p )',
}

fe = {
    'chunk_size' : 1000
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
        from sfepy.base.base import select_by_names
        from sfepy.fem import eval_term_op

        ok = True

        volumes = {}
        avg = 0.0
        for key, term in expressions.items():
            val = self.problem.evaluate(term)

            volumes[key] = val
            avg += val

        avg /= len(volumes)

        for key, val in volumes.items():
            err = nm.abs( avg - val ) / nm.abs( avg )
            _ok = err < 1e-12
            self.report( '"'"%s"'" - relative volume difference: %e -> %s'\
                         % (key, err, _ok) )
            ok = ok and _ok

        return ok

