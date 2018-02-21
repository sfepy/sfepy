from __future__ import absolute_import
import six
input_names = {'TL': '../examples/large_deformation/hyperelastic.py',
               'UL': '../examples/large_deformation/hyperelastic_ul.py',
               'ULM': '../examples/large_deformation/hyperelastic_ul_up.py'}
output_name_trunk = 'test_hyperelastic_'

from sfepy.base.testing import TestCommon
from tests_basic import NLSStatus

class Test(TestCommon):

    @staticmethod
    def from_conf(conf, options):
        return Test(conf = conf, options = options)

    def test_solution(self):

        from sfepy.base.base import Struct
        from sfepy.base.conf import ProblemConf, get_standard_keywords
        from sfepy.applications import solve_pde, assign_standard_hooks
        import numpy as nm
        import os.path as op

        solutions = {}
        ok = True

        for hp, pb_filename in six.iteritems(input_names):

            required, other = get_standard_keywords()
            input_name = op.join(op.dirname(__file__), pb_filename)
            test_conf = ProblemConf.from_file(input_name, required, other)

            name = output_name_trunk + hp
            solver_options = Struct(output_filename_trunk=name,
                                    output_format='vtk',
                                    save_ebc=False, save_ebc_nodes=False,
                                    save_regions=False,
                                    save_regions_as_groups=False,
                                    save_field_meshes=False,
                                    solve_not=False)
            assign_standard_hooks(self, test_conf.options.get, test_conf)

            self.report( 'hyperelastic formulation: %s' % (hp, ) )

            status = NLSStatus(conditions=[])

            pb, state = solve_pde(test_conf,
                                  solver_options,
                                  status=status,
                                  output_dir=self.options.out_dir,
                                  step_hook=self.step_hook,
                                  post_process_hook=self.post_process_hook,
                                  post_process_hook_final=self.post_process_hook_final)

            converged = status.nls_status.condition == 0
            ok = ok and converged

            solutions[hp] = state.get_parts()['u']
            self.report('%s solved' % input_name)

        rerr = 1.0e-3
        aerr = nm.linalg.norm(solutions['TL'], ord=None) * rerr

        self.report('allowed error: rel = %e, abs = %e' % (rerr, aerr))
        ok = ok and self.compare_vectors(solutions['TL'], solutions['UL'],
                                         label1='TLF',
                                         label2='ULF',
                                         allowed_error=rerr)

        ok = ok and self.compare_vectors(solutions['UL'], solutions['ULM'],
                                         label1='ULF',
                                         label2='ULF_mixed',
                                         allowed_error=rerr)

        return ok
