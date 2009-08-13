from sfepy.base.base import *
from sfepy.fem.mesh import find_refinement
from sfepy.fem.evaluate import eval_term
from sfepy.fem.fields import Fields
from sfepy.fem.variables import Variables
from sfepy.solvers import Solver

##
# 19.09.2007, c
def _make_full_and_save( filename, pb, svec, var_name ):
    from sfepy.base.ioutils import write_bb
    vec0 = pb.variables.make_full_vec( svec, var_name )
    vec = pb.variables[var_name].extend_data( vec0[:,nm.newaxis], pb.domain.n_nod )
    fd = open( filename, 'w' )
    write_bb( fd, vec, 2 )
    fd.close()


##
# 03.09.2007, c
class ProlongatorToQP( Struct ):

    ##
    # c: 03.09.2007, r: 13.02.2008
    def __init__( self, fpb, cpb, options, test_name, var_name, debug = False,
                  **kwargs ):
        Struct.__init__( self, **kwargs )

        self.debug = debug
        self.test_name, self.var_name = test_name, var_name

        eps = fpb.domain.get_diameter() * 1e-8
        coors = fpb.domain.get_mesh_coors()
        coarse_coors = cpb.domain.get_mesh_coors()

        tt = time.clock()
        aux = find_refinement( coarse_coors, cpb.domain.get_conns(),
                              coors, fpb.domain.get_conns(), eps )
        print time.clock() - tt
        self.emaps, self.iemaps, self.n_f_els = aux

        cfield = cpb.variables[var_name].field
        coarse_interps = cfield.get_interpolants()

        tt = time.clock()
        ffield = fpb.variables[var_name].field
        pbase = ffield.get_prolong_base( coors, cpb.domain, self.iemaps,
                                       coarse_interps, options.suppress_errors )
        print time.clock() - tt

        self.pbase = pbase

        region_name = ffield.bases[0][0] # !!!
        print region_name

#        fpb.time_update( conf_ebc = {}, conf_epbc = {} )
        dummy = fpb.create_state_vector()
        
        mtx_ff = fpb.variables.create_matrix_graph( var_names = [var_name] )
        print mtx_ff.shape

        mtx_term = 'dw_mass_scalar.%s( %s, %s )'\
                  % (region_name, test_name, var_name)
        mtx_ff = eval_term( dummy, mtx_term,
                          fpb.domain, fpb.regions, fpb.variables,
                          fpb.materials, chunk_size = fpb.domain.mesh.n_el,
                          dw_mode = 'matrix', tangent_matrix = mtx_ff )
##         import sfepy.base.plotutils as plu
##         plu.spy( mtx_ff, eps = 1e-8 )
##         plu.pylab.show()

        indx = fpb.variables.adi.indx[var_name]
        print indx
        self.mtx_ff = mtx_ff.get_submatrix( indx, indx )
        print self.mtx_ff.shape
##         import sfepy.base.plotutils as plu
##         plu.spy( self.mtx_ff, eps = 1e-8 )
##         plu.pylab.show()
        
##         vv = nla.svdvals( mtx_ff.toarray() )
##         smin, smax = vv.min(), vv.max()
##         cond = smax / smin
##         print smin, smax, cond

        vf = fpb.variables[test_name]
        vc = cpb.variables[var_name]
        fields = Fields.from_field_list( [vf.field, vc.field],
                                       fpb.fields.qp_coors,
                                       names = ['fine', 'coarse'] )

        vconf = {
            var_name  : ('field', 'unknown', 'coarse', vc.dof_types, 0),
            test_name : ('field', 'test', 'fine', vf.dof_types, var_name),
        }
        variables = Variables.from_conf( vconf, fields )
        variables.setup_dof_info( make_virtual = True )
        variables.setup_dof_conns( make_virtual = True )
        variables.equation_mapping( fpb.conf.ebc, fpb.conf.epbc,
                                   cpb.regions, None, fpb.conf.funmod,
                                   vregions = fpb.regions )
#        variables.equation_mapping( {}, {}, None, None, None )
        variables.setup_a_dof_conns()
        variables.fix_coarse_grid_a_dof_conns( self.iemaps, var_name )
##         print variables.adof_conns
##         print variables.avdof_conns
        mtx_fc = variables.create_matrix_graph()

        rhs_term = 'dw_mass_scalar_fine_coarse.%s( %s, %s, iemaps, pbase )'\
                  % (region_name, test_name, var_name)
#        debug()
        mtx_fc = eval_term( dummy, rhs_term,
                          fpb.domain, fpb.regions, variables,
                          fpb.materials, chunk_size = fpb.domain.mesh.n_el,
                          dw_mode = 'matrix', tangent_matrix = mtx_fc,
                          iemaps = self.iemaps, pbase = self.pbase )
        self.mtx_fc = mtx_fc
        print self.mtx_fc.shape

        assert_( self.mtx_ff.shape[0] == self.mtx_ff.shape[1] )
        assert_( self.mtx_ff.shape[0] == self.mtx_fc.shape[0] )

        self.ls_conf = fpb.get_solver_conf( 'ls' )
        if self.debug:
            self.cpb = cpb
            self.fpb = fpb

##             print variables.di
##             print variables.vdi
##             print variables.adi
##             print variables.avdi
##             pause()
##             import sfepy.base.plotutils as plu
##             plu.spy( mtx_fc, eps = 1e-8 )
##             plu.pylab.show()
##             debug()

#        fpb.time_update() # Restore EBC.

    ##
    # c: 05.09.2007, r: 13.02.2008
    def __call__( self, mtx_cb ):
        ls = Solver.any_from_conf( self.ls_conf, mtx = self.mtx_ff )
        
        assert_( self.mtx_fc.shape[1] == mtx_cb.shape[0] )

        mtx_fb = nm.empty( (self.mtx_ff.shape[0], mtx_cb.shape[1]),
                          dtype = mtx_cb.dtype )
        for ii in xrange( mtx_cb.shape[1] ):
            mtx_fb[:,ii] = ls( self.mtx_fc * mtx_cb[:,ii] )

##         import pylab as p
##         p.subplot( 211 )
##         p.plot( mtx_cb[:,0] )
##         p.subplot( 212 )
##         p.plot( mtx_fb[:,0] )
##         p.show()
        if self.debug:
            _make_full_and_save( self.cpb.domain.mesh.name + '.bb',
                              self.cpb, mtx_cb[:,0], self.var_name )
            _make_full_and_save( self.fpb.domain.mesh.name + '.bb',
                              self.fpb, mtx_fb[:,0], self.var_name )
            debug()
            
        return mtx_fb

