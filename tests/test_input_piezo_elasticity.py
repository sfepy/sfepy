from __future__ import absolute_import
input_name = '../examples/multi_physics/piezo_elasticity.py'
output_name = 'test_piezo_elasticity.vtk'


from tests_basic import TestInput

class Test( TestInput ):

    def from_conf( conf, options ):
        return TestInput.from_conf( conf, options, cls = Test )
    from_conf = staticmethod( from_conf )

    def test_ebc( self ):
        import numpy as nm
        from sfepy.discrete import Problem

        pb = Problem.from_conf(self.test_conf)
        pb.time_update()

        vvs = pb.get_variables()
        setv = vvs.set_state_part
        make_full = vvs.make_full_vec

        svec_u = nm.ones( (vvs.adi.n_dof['u'],), dtype = nm.float64 )
        svec_phi = nm.empty( (vvs.adi.n_dof['phi'],), dtype = nm.float64 )
        svec_phi.fill( 2.0 )

        svec = vvs.create_stripped_state_vector()
        setv( svec, svec_u, 'u', stripped = True )
        setv( svec, svec_phi, 'phi', stripped = True )

        vec = make_full( svec )

        ii_u = vvs.di.indx['u'].start + vvs['u'].eq_map.eqi
        ii_phi = vvs.di.indx['phi'].start + vvs['phi'].eq_map.eqi

        ok_ebc = vvs.has_ebc( vec )
        ok_u = nm.all( vec[ii_u] == svec_u )
        ok_phi = nm.all( vec[ii_phi] == svec_phi )

        msg = '%s: %s'
        self.report( msg % ('ebc', ok_ebc) )
        self.report( msg % ('u', ok_u) )
        self.report( msg % ('phi', ok_phi) )

        ok = ok_ebc and ok_u and ok_phi

        return ok
