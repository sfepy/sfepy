import os.path as op

from sfepy import data_dir
from sfepy.fem.periodic import match_y_line

filename_mesh = data_dir + '/meshes/2d/square_unit_tri.mesh'

materials = {
    'm' : ({'K' : [[3.0, 0.1], [0.3, 1.0]]},),
}

fields = {
    'pressure' : ('real', 'scalar', 'Omega', 2),
}

variables = {
    'p' : ('unknown field', 'pressure', 0),
    'q' : ('test field', 'pressure', 'p'),
    'r' : ('parameter field', 'pressure', 'p'),
}

regions = {
    'Omega' : ('all', {}),
    'Left' : ('nodes in (x < -0.499)', {}),
    'LeftStrip' : ('nodes in (x < -0.499) & (y > -0.199) & (y < 0.199)', {}),
    'LeftFix' : ('r.Left -n r.LeftStrip', {}),
    'Right' : ('nodes in (x > 0.499)', {}),
    'RightStrip' : ('nodes in (x > 0.499) & (y > -0.199) & (y < 0.199)', {}),
    'RightFix' : ('r.Right -n r.RightStrip', {}),
}

ebcs = {
    't_left' : ('LeftFix', {'p.0' : 5.0}),
    't_right' : ('RightFix', {'p.0' : 0.0}),
}

epbcs = {
    'periodic_x' : (['LeftStrip', 'RightStrip'],
                    {'p.0' : 'p.0'}, 'match_y_line'),
}

functions = {
    'match_y_line' : (match_y_line,),
}

equations = {
    'eq' : """dw_diffusion.2.Omega( m.K, q, p ) = 0"""
}

solvers = {
    'ls' : ('ls.scipy_direct', {}),
    'newton' : ('nls.newton',
                {'i_max'      : 1,
                 'eps_a'      : 1e-10,
    }),
}

from sfepy.base.testing import TestCommon

class Test(TestCommon):

    @staticmethod
    def from_conf(conf, options):
        from sfepy.applications import solve_pde

        problem, state = solve_pde(conf)
        name = op.join(options.out_dir,
                       op.splitext(op.basename(__file__))[0] + '.vtk')
        problem.save_state(name, state)

        test = Test(problem=problem, state=state, conf=conf, options=options)
        return test

    def test_save_ebc(self):
        name = op.join(self.options.out_dir,
                       op.splitext(op.basename(__file__))[0])
        self.problem.save_ebc(name + '_ebcs_f.vtk', force=True)
        self.problem.save_ebc(name + '_ebcs.vtk', default=-1, force=False)

        return True
