# c: 06.02.2008, r: 06.02.2008
input_name = 'input/time_poisson.py'
output_name_trunk = 'test_time_poisson'

from testsBasic import TestInput

##
# c: 06.02.2008
class Test( TestInput ):

    ##
    # c: 06.02.2008, r: 06.02.2008
    def from_conf( conf, options ):
        return TestInput.from_conf( conf, options, cls = Test )
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
        options = Struct( output_file_name_trunk = name,
                          save_ebc = False, save_regions = False,
                          save_field_meshes = False, save_region_field_meshes = False,
                          solve_not = False )
        dpb, vec_dp, data = solve_direct( self.test_conf, options )

        self.report( '%s solved' % self.conf.input_name )

        return True
