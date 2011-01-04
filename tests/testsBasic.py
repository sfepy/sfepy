"""
This module is not a test file. It contains classes grouping some common
functionality, that is used in several test files.
"""
from sfepy.base.base import IndexedStruct
from sfepy.base.testing import TestCommon
import os.path as op

class NLSStatus(IndexedStruct):
    """
    Custom nonlinear solver status storing stopping condition of all
    time steps.
    """
    def __setitem__(self, key, val):
        IndexedStruct.__setitem__(self, key, val)
        if key == 'condition':
            self.conditions.append(val)

class TestInput(TestCommon):
    """Test that an input file works. See test_input_*.py files."""

    @staticmethod
    def from_conf(conf, options, cls=None):
        from sfepy.base.base import Struct
        from sfepy.base.conf import ProblemConf, get_standard_keywords
        from sfepy.applications.simple_app import assign_standard_hooks

        required, other = get_standard_keywords()
        input_name = op.join(op.dirname(__file__), conf.input_name)
        test_conf = ProblemConf.from_file(input_name, required, other)

        if cls is None:
            cls = TestInput
        test = cls(test_conf=test_conf, conf=conf, options=options)

        assign_standard_hooks(test, test_conf.options.get_default_attr,
                              test_conf.funmod)

        name = op.join(test.options.out_dir, test.get_output_name_trunk())
        test.solver_options = Struct(output_filename_trunk = name,
                                     output_format ='vtk',
                                     save_ebc = False, save_regions = False,
                                     save_regions_as_groups = False,
                                     save_field_meshes = False,
                                     solve_not = False)

        return test

    def get_output_name_trunk(self):
        return op.splitext(op.split(self.conf.output_name)[1])[0]

    def check_conditions(self, conditions):
        ok = (conditions == 0).all()
        if not ok:
            self.report('nls stopping conditions:')
            self.report(conditions)
        return ok

    def test_input(self):
        import numpy as nm
        from sfepy.solvers.generic import solve_direct

        self.report('solving %s...' % self.conf.input_name)

        status = NLSStatus(conditions=[])

        out = solve_direct(self.test_conf,
                           self.solver_options,
                           step_hook=self.step_hook,
                           post_process_hook=self.post_process_hook,
                           post_process_hook_final=self.post_process_hook_final,
                           nls_status=status)
        self.report('%s solved' % self.conf.input_name)

        ok = self.check_conditions(nm.array(status.conditions))

        return ok

class TestInputEvolutionary(TestInput):

    @staticmethod
    def from_conf(conf, options, cls=None):
        if cls is None:
            cls = TestInputEvolutionary

        return TestInput.from_conf(conf, options, cls=cls)

    def get_output_name_trunk(self):
        return self.conf.output_name_trunk

##
# 03.10.2007, c
class TestLCBC( TestCommon ):
    """Test linear combination BC. See test_lcbc_*.py files."""

    ##
    # 03.10.2007, c
    def from_conf( conf, options ):
        return TestLCBC( conf = conf, options = options )
    from_conf = staticmethod( from_conf )

    ##
    # c: 03.10.2007, r: 08.02.2008
    def test_linear_rigid_body_bc( self ):
        import scipy
        if scipy.version.version == "0.6.0":
            # This test uses a functionality implemented in scipy svn, which is
            # missing in scipy 0.6.0
            return True
        from sfepy.base.base import Struct
        from sfepy.solvers.generic import solve_stationary
        from sfepy.base.base import IndexedStruct

        status = IndexedStruct()
        problem, state = solve_stationary(self.conf, nls_status=status)
        ok = status.condition == 0
        self.report( 'converged: %s' % ok )
        out = state.create_output_dict()

        strain = problem.evaluate('de_cauchy_strain.i1.Y( u )', mode='el_avg')
        out['strain'] = Struct(name='output_data',
                               mode='cell', data=strain, dofs=None)

        name = op.join( self.options.out_dir,
                        op.split( self.conf.output_name )[1] )
        problem.domain.mesh.write( name, io = 'auto', out = out )

        ##
        # Check if rigid body displacements are really rigid should go here.

        return ok
