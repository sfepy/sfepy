import pytest

import sfepy.base.testing as tst

input_names = {'TL': 'examples/large_deformation/hyperelastic.py',
               'UL': 'examples/large_deformation/hyperelastic_ul.py',
               'ULM': 'examples/large_deformation/hyperelastic_ul_up.py'}
output_name_trunk = 'test_hyperelastic_'

@pytest.mark.slow
def test_solution(output_dir):
    import sfepy
    from sfepy.base.base import IndexedStruct, Struct
    from sfepy.applications import solve_pde
    import numpy as nm
    import os.path as op

    solutions = {}
    ok = True

    for hp, pb_filename in input_names.items():

        input_name = op.join(sfepy.base_dir, pb_filename)

        name = output_name_trunk + hp
        solver_options = Struct(output_filename_trunk=name,
                                output_format='vtk',
                                save_ebc=False, save_ebc_nodes=False,
                                save_regions=False,
                                save_regions_as_groups=False,
                                solve_not=False)

        tst.report('hyperelastic formulation: %s' % hp)

        status = IndexedStruct(nls_status=tst.NLSStatus(conditions=[]))

        pb, state = solve_pde(input_name, solver_options, status=status,
                              output_dir=output_dir)
        converged = status.nls_status.condition == 0
        ok = ok and converged

        solutions[hp] = state.get_state_parts()['u']
        tst.report('%s solved' % input_name)

    rerr = 1.0e-3
    aerr = nm.linalg.norm(solutions['TL'], ord=None) * rerr

    tst.report('allowed error: rel = %e, abs = %e' % (rerr, aerr))
    ok = ok and tst.compare_vectors(solutions['TL'], solutions['UL'],
                                    label1='TLF',
                                    label2='ULF',
                                    allowed_error=rerr)

    ok = ok and tst.compare_vectors(solutions['UL'], solutions['ULM'],
                                    label1='ULF',
                                    label2='ULF_mixed',
                                    allowed_error=rerr)

    assert ok
