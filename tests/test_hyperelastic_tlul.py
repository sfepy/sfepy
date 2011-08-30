input_names = {'TL': '../examples/large_deformation/hyperelastic.py',
               'UL': '../examples/large_deformation/hyperelastic_ul.py'}
output_name_trunk = 'test_hyperelastic_'

from sfepy.base.testing import TestCommon
from testsBasic import NLSStatus

class Test(TestCommon):

    @staticmethod
    def from_conf(conf, options):
        return Test(conf = conf, options = options)

    def test_solution(self):

        from sfepy.base.base import Struct
        from sfepy.base.conf import ProblemConf, get_standard_keywords
        from sfepy.applications.simple_app import assign_standard_hooks
        from sfepy.solvers.generic import solve_direct
        import numpy as nm
        import os.path as op

        solutions = {}
        ok = True

        for hp, pb_filename in input_names.iteritems():

            required, other = get_standard_keywords()
            input_name = op.join(op.dirname(__file__), pb_filename)
            test_conf = ProblemConf.from_file(input_name, required, other)

            solver_options = Struct(output_filename_trunk = output_name_trunk + hp,
                                    output_format ='vtk',
                                    save_ebc = False, save_regions = False,
                                    save_regions_as_groups = False,
                                    save_field_meshes = False,
                                    solve_not = False)
            assign_standard_hooks(self, test_conf.options.get_default_attr,
                                  test_conf)

            self.report( 'hyperelastic formulation: %s' % (hp, ) )

            status = NLSStatus(conditions=[])

            pb, state = solve_direct(test_conf,
                                     solver_options,
                                     step_hook=self.step_hook,
                                     post_process_hook=self.post_process_hook,
                                     post_process_hook_final=self.post_process_hook_final,
                                     nls_status=status)

            converged = status.condition == 0
            ok = ok and converged

            solutions[hp] = state.get_parts()['u']
            self.report('%s solved' % input_name)

        rerr = 1.0e-2;
        aerr = nm.linalg.norm(solutions['TL'], ord=None) * rerr

        self.report('allowederror: rel = %e, abs = %e' % (rerr, aerr) )
        ok = ok and self.compare_vectors( solutions['TL'], solutions['UL'],
                                          label1 = 'TLF',
                                          label2 = 'ULF',
                                          allowed_error = rerr)

        return ok
