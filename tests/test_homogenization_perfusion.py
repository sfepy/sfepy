from __future__ import print_function
from __future__ import absolute_import
import six
input_name = '../examples/homogenization/perfusion_micro.py'

from sfepy.base.testing import TestCommon

class Test(TestCommon):

    @staticmethod
    def from_conf(conf, options):
        return Test(conf = conf, options = options)

    def compare_scalars(s1, s2, l1= 's1', l2 = 's2',
                        allowed_error = 1e-8):

        diff  = abs(s1 - s2)
        TestCommon.report( '|%s - %s|: %e' % (l1, l2, diff))
        if diff > allowed_error:
            return False
        else:
            return True
    compare_scalars = staticmethod(compare_scalars)

    def test_solution(self):

        from sfepy.base.base import Struct
        from sfepy.base.conf import ProblemConf, get_standard_keywords
        from sfepy.homogenization.homogen_app import HomogenizationApp
        #import numpy as nm
        import os.path as op

        ok = True

        required, other = get_standard_keywords()
        required.remove('equations')
        print(input_name)
        full_name = op.join(op.dirname(__file__), input_name)
        test_conf = ProblemConf.from_file(full_name, required, other)

        options = Struct(output_filename_trunk=None,
                         save_ebc=False,
                         save_ebc_nodes=False,
                         save_regions=False,
                         save_field_meshes=False,
                         save_regions_as_groups=False,
                         solve_not=False)

        test_conf.options['output_dir'] = './output-tests'
        app = HomogenizationApp(test_conf, options, 'homogen:' )
        coefs = app()

        aerr = 1.0e-9
        self.report('allowed error: abs = %e' % (aerr, ))

        # G^A = G^B ?
        ok = ok and self.compare_scalars(coefs.GA, coefs.GB,\
                                         'G^A', 'G^B', aerr)

        # F^{A+} + F^{B+} = -1/h \int_{\partial_+Y_m} ?
        aux = 1.0 / test_conf.param_h * coefs.Volume_bYMp
        ok = ok and self.compare_scalars(coefs.FpA + coefs.FpB, -aux,
                                         'F^{A+} + F^{B+}', '-bYM^+', aerr)

        # F^{A-} + F^{B-} = -1/h \int_{\partial_-Y_m} ?
        aux = 1.0 / test_conf.param_h * coefs.Volume_bYMm
        ok = ok and self.compare_scalars(coefs.FmA + coefs.FmB, -aux,
                                         'F^{A-} + F^{B-}', '-bYM^-', aerr)

        # symmetry of H ?
        ok = ok and self.compare_scalars(coefs.Hpm, coefs.Hmp,
                                         'H^{+-}', 'H^{-+}', aerr)

        # E = -F ?
        ok = ok and self.compare_scalars(coefs.EmA, -coefs.FmA,
                                         'E^{A-}', '-F^{A-}',aerr)
        ok = ok and self.compare_scalars(coefs.EpA, -coefs.FpA,
                                         'E^{A+}', '-F^{A+}',aerr)
        ok = ok and self.compare_scalars(coefs.EmB, -coefs.FmB,
                                         'E^{B-}', '-F^{B-}',aerr)
        ok = ok and self.compare_scalars(coefs.EpB, -coefs.FpB,
                                         'E^{B+}', '-F^{B+}',aerr)

        # S = S_test ?
        coefsd = coefs.to_dict()
        compare = []
        for ii in six.iterkeys(coefsd):
            if 'S_test' in ii:
                ch = ii[6]
                io = ii[-1]
                compare.append((ii, 'S%s_%s' % (ch, io)))

        for s1, s2 in compare:
            ok = ok and self.compare_vectors(coefsd[s1], -coefsd[s2],
                                             label1='S_test', label2='S',
                                             allowed_error=aerr)

        return ok
