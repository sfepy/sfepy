import os.path as op

from sfepy import data_dir

from sfepy.base.testing import TestCommon

class Test(TestCommon):

    @staticmethod
    def from_conf(conf, options):
        from sfepy.applications import solve_pde

        problem, state = solve_pde(conf, output_dir=options.out_dir)

        test = Test(problem=problem, state=state, conf=conf, options=options)
        return test

    def test_save_ebc(self):
        name = op.join(self.options.out_dir,
                       op.splitext(op.basename(__file__))[0])
        self.problem.save_ebc(name + '_ebcs_f.vtk', force=True)
        self.problem.save_ebc(name + '_ebcs.vtk', default=-1, force=False)

        return True
