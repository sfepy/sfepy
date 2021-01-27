"""
Test computing volumes by volume or surface integrals.
"""
from __future__ import absolute_import
from sfepy import data_dir

filename_mesh = data_dir + '/meshes/3d/elbow.mesh'

fields = {
    'scalar' : ('real', 'scalar', 'Omega', 1),
    'vector' : ('real', 'vector', 'Omega', 1),
}

integrals = {
    'i' : 2,
}

regions = {
    'Omega' : 'all',
    'Gamma' : ('vertices of surface', 'facet'),
}

expressions = {
    'volume_p' : 'd_volume.i.Omega(p)',
    'volume_u' : 'd_volume.i.Omega(u)',
    'surface_p' : 'd_volume_surface.i.Gamma(p)',
    'surface_u' : 'd_volume_surface.i.Gamma(u)',
}

import numpy as nm
from sfepy.base.testing import TestCommon

class Test(TestCommon):

    @staticmethod
    def from_conf(conf, options):
        from sfepy.discrete import Problem

        problem = Problem.from_conf(conf, init_equations=False)
        test = Test(problem=problem, conf=conf, options=options)
        return test

    def test_volume(self):
        from sfepy.discrete import FieldVariable

        ok = True

        field_map = {'u' : 'vector', 'p' : 'scalar'}

        volumes = {}
        avg = 0.0
        for key, term in expressions.items():
            var_name = key[-1]
            field = self.problem.fields[field_map[var_name]]
            var = FieldVariable(var_name, 'parameter', field,
                                primary_var_name='(set-to-None)')

            val = self.problem.evaluate(term, **{var_name : var})

            volumes[key] = val
            avg += val

        avg /= len(volumes)

        for key, val in volumes.items():
            err = nm.abs(avg - val) / nm.abs(avg)
            _ok = err < 1e-12
            self.report('"'"%s"'" - volume: %e' % (key, val))
            self.report('"'"%s"'" - relative volume difference: %e -> %s'
                        % (key, err, _ok))
            ok = ok and _ok

        return ok

    def test_volume_tl(self):
        from sfepy.discrete import FieldVariable

        fu = self.problem.fields['vector']
        fq = self.problem.fields['scalar']

        var_u = FieldVariable('u', 'parameter', fu,
                              primary_var_name='(set-to-None)')
        var_q = FieldVariable('q', 'test', fq,
                              primary_var_name='(set-to-None)')

        var_u.set_data(nm.linspace(0, 0.004, var_u.n_dof))

        vval = self.problem.evaluate('dw_tl_volume.i.Omega( q, u )',
                                     term_mode='volume', q=var_q, u=var_u)

        sval = self.problem.evaluate('d_tl_volume_surface.i.Gamma( u )',
                                     u=var_u)

        ok = abs(vval - sval) < 1e-14

        self.report('TL: by volume: %e == by surface: %e -> %s' %
                    (vval, sval, ok))

        return ok
