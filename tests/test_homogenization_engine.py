from __future__ import absolute_import
from sfepy.base.testing import TestCommon

class Test(TestCommon):

    @staticmethod
    def from_conf(conf, options):
        test = Test(conf=conf, options=options)
        return test

    def test_dependencies(self):
        from sfepy.homogenization.engine import HomogenizationEngine

        get_deps = HomogenizationEngine.get_sorted_dependencies

        coefs = {'A' : {'requires' : ['a', 'd', 'c.B']},
                 'B' : {'requires' : ['b']}}
        requirements = {'a' : {'requires' : ['b', 'c']},
                        'b' : {'requires' : ['c']},
                        'c' : {},
                        'd' : {'requires' : ['b', 'a']}}

        deps = get_deps(requirements, coefs, None)
        ok = ((deps == ['c', 'b', 'a', 'd', 'c.B', 'c.A'])
              or (deps == ['c', 'b', 'c.B', 'a', 'd', 'c.A'])
              or (deps == ['c', 'b', 'a', 'c.B', 'd', 'c.A']))
        self.report(deps, ':', ok)

        coefs['B']['requires'] = ['b', 'c.A']

        try:
            deps = get_deps(requirements, coefs, None)

        except ValueError as err:
            self.report('detected:', str(err))
            _ok = 'circular requirement "c.' in str(err)

        else:
            _ok = False

        self.report('circular dependency detection 1:', _ok)
        ok = ok and _ok

        coefs['B']['requires'] = ['b']
        requirements['c']['requires'] = ['d']

        try:
            deps = get_deps(requirements, coefs, None)

        except ValueError as err:
            self.report('detected:', str(err))
            _ok = 'circular requirement' in str(err)

        else:
            _ok = False

        self.report('circular dependency detection 2:', _ok)
        ok = ok and _ok

        return ok
