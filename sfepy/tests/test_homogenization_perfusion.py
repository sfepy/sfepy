import pytest

import sfepy.base.testing as tst

input_name = 'examples/homogenization/perfusion_micro.py'

def compare_scalars(s1, s2, l1= 's1', l2 = 's2',
                    allowed_error = 1e-8):

    diff  = abs(s1 - s2)
    tst.report( '|%s - %s|: %e' % (l1, l2, diff))
    if diff > allowed_error:
        return False
    else:
        return True

@pytest.mark.slow
def test_solution(output_dir):
    import sfepy
    from sfepy.base.base import Struct
    from sfepy.base.conf import ProblemConf, get_standard_keywords
    from sfepy.homogenization.homogen_app import HomogenizationApp
    import os.path as op

    ok = True

    required, other = get_standard_keywords()
    required.remove('equations')
    full_name = op.join(sfepy.base_dir, input_name)
    test_conf = ProblemConf.from_file(full_name, required, other)

    options = Struct(output_filename_trunk=None,
                     save_ebc=False,
                     save_ebc_nodes=False,
                     save_regions=False,
                     save_regions_as_groups=False,
                     solve_not=False)

    test_conf.options['output_dir'] = output_dir
    app = HomogenizationApp(test_conf, options, 'homogen:' )
    coefs = app()

    aerr = 1.0e-9
    tst.report('allowed error: abs = %e' % aerr)

    # G^A = G^B ?
    ok = ok and compare_scalars(coefs.GA, coefs.GB,\
                                'G^A', 'G^B', aerr)

    # F^{A+} + F^{B+} = -1/h \int_{\partial_+Y_m} ?
    aux = 1.0 / test_conf.param_h * coefs.Volume_bYMp
    ok = ok and compare_scalars(coefs.FpA + coefs.FpB, -aux,
                                'F^{A+} + F^{B+}', '-bYM^+', aerr)

    # F^{A-} + F^{B-} = -1/h \int_{\partial_-Y_m} ?
    aux = 1.0 / test_conf.param_h * coefs.Volume_bYMm
    ok = ok and compare_scalars(coefs.FmA + coefs.FmB, -aux,
                                'F^{A-} + F^{B-}', '-bYM^-', aerr)

    # symmetry of H ?
    ok = ok and compare_scalars(coefs.Hpm, coefs.Hmp,
                                'H^{+-}', 'H^{-+}', aerr)

    # E = -F ?
    ok = ok and compare_scalars(coefs.EmA, -coefs.FmA,
                                'E^{A-}', '-F^{A-}',aerr)
    ok = ok and compare_scalars(coefs.EpA, -coefs.FpA,
                                'E^{A+}', '-F^{A+}',aerr)
    ok = ok and compare_scalars(coefs.EmB, -coefs.FmB,
                                'E^{B-}', '-F^{B-}',aerr)
    ok = ok and compare_scalars(coefs.EpB, -coefs.FpB,
                                'E^{B+}', '-F^{B+}',aerr)

    # S = S_test ?
    coefsd = coefs.to_dict()
    compare = []
    for ii in coefsd.keys():
        if 'S_test' in ii:
            ch = ii[6]
            io = ii[-1]
            compare.append((ii, 'S%s_%s' % (ch, io)))

    for s1, s2 in compare:
        ok = ok and tst.compare_vectors(coefsd[s1], -coefsd[s2],
                                        label1='S_test', label2='S',
                                        allowed_error=aerr)

    assert ok
