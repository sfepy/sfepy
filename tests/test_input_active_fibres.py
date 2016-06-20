from __future__ import absolute_import
input_name = '../examples/large_deformation/active_fibres.py'
output_name_trunk = 'test_active_fibres'

from tests_basic import TestInputEvolutionary

class Test(TestInputEvolutionary):

    @staticmethod
    def from_conf(conf, options):
        return TestInputEvolutionary.from_conf(conf, options, cls=Test)

    def check_conditions(self, conditions):
        """
        Special-case the first iteration, as the solver converges slowly there.
        """
        ok = (conditions[1:] == 0).all()
        ok = ok and (conditions[0] == 1)

        if not ok:
            self.report('nls stopping conditions:')
            self.report(conditions)

        return ok
