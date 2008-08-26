"""
This module is not a test file. It contains classes grouping some common
functionality, that is used in several test files.
"""

from sfepy.base.testing import TestCommon
import os.path as op

##
# 05.06.2007, c
class TestInput( TestCommon ):
    """Test that an input file works. See test_input_*.py files."""

    ##
    # c: 05.06.2007, r: 19.02.2008
    def from_conf( conf, options, cls = None ):
        from sfepy.base.conf import ProblemConf, get_standard_keywords

        required, other = get_standard_keywords()
        test_conf = ProblemConf.from_file( conf.input_name, required, other )

        if cls is None:
            cls = TestInput
        test = cls( test_conf = test_conf, conf = conf, options = options )

        return test
    from_conf = staticmethod( from_conf )

    ##
    # c: 05.06.2007, r: 08.02.2008
    def test_input( self ):
        from sfepy.solvers.generic import solve_stationary
        from sfepy.base.base import IndexedStruct

        self.report( 'solving %s...' % self.conf.input_name )
        status = IndexedStruct()
        dpb, vec_dp, data = solve_stationary( self.test_conf, nls_status = status )
        ok = status.condition == 0
        out = dpb.state_to_output( vec_dp )

        name = op.join( self.options.out_dir,
                        op.split( self.conf.output_name )[1] )
        dpb.save_state( name, vec_dp )
        self.report( '%s solved' % self.conf.input_name )

        return ok

class TestInputEvolutionary( TestInput ):

    ##
    # c: 06.02.2008, r: 06.02.2008
    def from_conf( conf, options ):
        return TestInput.from_conf( conf, options, cls = TestInputEvolutionary )
    from_conf = staticmethod( from_conf )

    ##
    # c: 06.02.2008, r: 06.02.2008
    def test_input( self ):
        """Does not verify anything!!!"""
        import os.path as op
        from sfepy.solvers.generic import solve_direct
        from sfepy.base.base import Struct

        self.report( 'solving %s...' % self.conf.input_name )
        name = op.join( self.options.out_dir, self.conf.output_name_trunk )
        options = Struct( output_filename_trunk = name,
                          save_ebc = False, save_regions = False,
                          save_field_meshes = False,
                          save_region_field_meshes = False,
                          solve_not = False )
        dpb, vec_dp, data = solve_direct( self.test_conf, options )

        self.report( '%s solved' % self.conf.input_name )

        return True
    

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
        from sfepy.fem.evaluate import eval_term_op

        status = IndexedStruct()
        problem, vec, data = solve_stationary( self.conf,
                                              nls_status = status )
        ok = status.condition == 0
        self.report( 'converged: %s' % ok )
        out = problem.state_to_output( vec )

        strain = eval_term_op( vec, 'de_cauchy_strain.i1.Y( u )', problem )
        out['strain'] = Struct( name = 'output_data',
                                mode = 'cell', data = strain,
                                dof_types = None )

        name = op.join( self.options.out_dir,
                        op.split( self.conf.output_name )[1] )
        problem.domain.mesh.write( name, io = 'auto', out = out )

        ##
        # Check if rigid body displacements are really rigid should go here.

        return ok
