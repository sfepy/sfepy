import os.path as op

from sfepy.base.base import *
from sfepy.base.progressbar import MyBar
import sfepy.fem.evaluate as eva
import freeFormDef as ffd
from sfepy.terms import CharacteristicFunction

##
# 27.07.2006, c
# 01.11.2007
def set_state_to_vars( variables, var_map, vec_dp, vec_ap = None ):
    # Set admissible state variables.
    for dvar, avar in var_map.iteritems():
        print '>>>>>>>>>>>>>>', dvar, avar
        variables.non_state_data_from_state( dvar, vec_dp, avar )
    if vec_ap is not None:
        # Set adjoint state variables.
        variables.data_from_state( vec_ap )

##
# c: 15.10.2007, r: 15.04.2008
def update_mesh( shape_opt, pb, design ):
    output( 'mesh update (%d)!' % shape_opt.cache.i_mesh )
    ##
    # Update mesh.
    shape_opt.dsg_vars.val = design
    shape_opt.sp_boxes.set_control_points( shape_opt.dsg_vars )
    coors = shape_opt.sp_boxes.interp_coordinates()
    # Do not update state, so that even warped mesh gets saved...
    pb.set_mesh_coors( coors, update_state = False )

    pb.domain.mesh.write( op.join( shape_opt.save_dir, 'design.%03d.mesh' )\
                          % shape_opt.cache.i_mesh, io = 'auto' )
    shape_opt.cache.i_mesh += 1
    try:
        pb.set_mesh_coors( coors, update_state = True )
    except:
        output( '...failed!' )
        return False
    return True

##
# c: 19.04.2006, r: 15.04.2008
def solve_problem_for_design(problem, design, shape_opt, opts,
                             use_cache=True, is_mesh_update=True):
    """use_cache == True means direct problem..."""
    pb = problem
    
    if use_cache and nm.allclose( design, shape_opt.cache.design,
                                 rtol = opts.eps_r, atol = opts.eps_a ):
        output( 'cache!' )
        vec = shape_opt.cache.vec
    else:
        if is_mesh_update:
            ok = update_mesh( shape_opt, pb, design )
            if not ok: return None
            
        # Solve direct/adjoint problem.
        vec = problem.solve()

        if use_cache:
            shape_opt.cache.vec = vec.copy()
            shape_opt.cache.design = design.copy()

            if shape_opt.save_iter_sols:
                pb.save_state( op.join( shape_opt.save_dir, 'direct.%03d.vtk' )\
                              % (shape_opt.cache.i_mesh - 1), vec )
                
            if shape_opt.save_control_points:
                aux = shape_opt.sp_boxes.create_mesh_from_control_points()
                aux.write( op.join( shape_opt.save_dir, 'cp.%03d.vtk' )\
                           % (shape_opt.cache.i_mesh - 1), io = 'auto' )

            if shape_opt.save_dsg_vars:
                filename = op.join( shape_opt.save_dir, 'dsg.%03d.txt' )\
                           % (shape_opt.cache.i_mesh - 1)
                shape_opt.dsg_vars.val.tofile( filename, ' ' )
##     print use_cache, is_mesh_update
##     print vec
    return vec

##
# created:       19.04.2006
# last revision: 21.12.2007
def obj_fun( design, shape_opt, dpb, apb, opts ):

#    print 'OF'
    data = {}
    vec_dp = solve_problem_for_design(dpb, design, shape_opt, opts )
    if vec_dp is None:
        return None
        
    val = shape_opt.obj_fun( vec_dp, dpb.conf, dpb.domain, dpb.variables,
                           dpb.materials, data = data )
#    print 'OF ok'
    return val

##
# created:       19.04.2006
# last revision: 21.12.2007
def obj_fun_grad( design, shape_opt, dpb, apb, opts ):

#    print 'OFG'
    data = {}
    vec_dp = solve_problem_for_design(dpb, design, shape_opt, opts )
    set_state_to_vars( apb.variables, shape_opt.var_map, vec_dp )
    vec_ap = solve_problem_for_design(apb, design, shape_opt, opts,
                                      use_cache=False, is_mesh_update=False)
##     print apb.materials['stabil']
##     pause()
    vec_sa = shape_opt.sensitivity( vec_dp, vec_ap, apb.conf, apb.domain,
                                  apb.variables, apb.materials, data = data )
#    print 'OFG ok'
    return vec_sa

##
# 13.04.2006, c
# 26.07.2006
# 20.03.2007
# 22.10.2007
# 24.10.2007
# 25.10.2007
def test_terms( idsgs, delta, shape_opt, vec_dp, vec_ap, pb ):

    dd = {}
    ccs = shape_opt.check_custom_sensitivity

##     ccs( 'd_sd_st_grad_div.Omega_D( stabil.gamma, w, u, Nu, mode )',
##          idsgs, delta, vec_dp, vec_ap, pb, dd )

##     ccs( 'd_sd_st_supg_c.Omega_D( stabil.delta, w, w, u, Nu, mode )',
##          idsgs, delta, vec_dp, vec_ap, pb, dd )

##     ccs( 'd_sd_st_pspg_c.Omega_D( stabil.tau, r, w, u, Nu, mode )',
##          idsgs, delta, vec_dp, vec_ap, pb, dd )

##     pb.variables.non_state_data_from_state( 'p', vec_dp, 'r' )

##     indx = pb.variables.get_indx( 'r' )
##     dummy = pb.create_state_vector()
##     mtx_pp = eva.eval_term_op( dummy, 'dw_st_pspg_p.Omega_D( stabil.tau, q, r )',
##                             pb, dw_mode = 'matrix' )[indx,indx]
##     v_p = vec_dp[indx]
##     v_r = vec_ap[indx]
##     pp1 = sc.dot( v_r.T * mtx_pp, v_p )
##     pp2 = eva.eval_term_op( vec_ap,
##                           'd_sd_st_pspg_p.Omega_D( stabil.tau, r, p, r, mode )',
##                           pb, mode = 0 )
##     print pp1, pp2
##     pause()
##     import sfepy.base.plotutils as plu
##     print mtx_pp
##     plu.spy( mtx_pp )
##     plu.pylab.show()
##     pause()
#    pb.time_update( None, conf_ebc = {} )
#    vec_dp.fill( 1.0 )
#    iop = eva.eval_term_op( vec_dp, 'dw_volume_integrate.Omega_C( q )',
#                          pb, dw_mode = 'vector' )
#    print iop, iop.shape
#    print sc.dot( vec_dp, iop )
#    oo = eva.eval_term_op( vec_dp, 'd_volume_integrate.Omega_C( p )', pb )
#    print oo
#    oo = eva.eval_term_op( vec_dp, 'd_volume.Omega_C( p )', pb )
#    print oo
#    debug()

    ccs( 'd_sd_st_pspg_p.i1.Omega( stabil.tau, r, p, Nu, mode )',
         idsgs, delta, vec_dp, vec_ap, pb, dd )
    ccs( 'd_sd_st_pspg_p.i1.Omega( stabil.tau, r, r, Nu, mode )',
         idsgs, delta, vec_dp, vec_ap, pb, dd )
    ccs( 'd_sd_st_pspg_p.i1.Omega( stabil.tau, p, p, Nu, mode )',
         idsgs, delta, vec_dp, vec_ap, pb, dd )

##     ccs( 'd_sd_test_pq.Omega_D( p, r, Nu, mode )',
##          idsgs, delta, vec_dp, vec_ap, pb, dd )

##     ccs( 'd_sd_div.Omega_D( u, r, Nu, mode )',
##          idsgs, delta, vec_dp, vec_ap, pb, dd )

##     ccs( 'd_sd_div.Omega_D( w, p, Nu, mode )',
##          idsgs, delta, vec_dp, vec_ap, pb, dd )

##     ccs( 'd_sd_div_grad.Omega_D( one, fluid, u, w, Nu, mode )',
##          idsgs, delta, vec_dp, vec_ap, pb, dd )

##     ccs( 'd_sd_convect.Omega_D( u, w, Nu, mode )',
##          idsgs, delta, vec_dp, vec_ap, pb, dd )

##
# 25.01.2006, c
class ShapeOptimFlowCase( Struct ):

    ##
    # c: 25.01.2006, r: 15.04.2008
    def from_conf( conf, domain ):
#        print conf_design

        opts = conf.options
        regions = domain.regions
        if opts.use_mesh_velocity:
            import tables as pt
            fd = pt.open_file( opts.mesh_velocity_filename, mode = 'r' )
            aux = fd.get_node( '/u' ).read()
            nu = nm.asarray( aux, dtype = nm.float64 )
            fd.close()

        else:
            nu = None

        sp_boxes = ffd.read_spline_box_hdf5( opts.ffd_spline_data )
        dsg_vars = ffd.read_dsg_vars_hdf5( opts.ffd_spline_data )
        dsg_vars.renumber_by_boxes( sp_boxes )
        dsg_vars.normalize_null_space_base()
        print dsg_vars.indx.shape
        print dsg_vars.null_space_b.shape

        control_region = regions[opts.control_domain]
        design_region = regions[opts.design_domain]

        from sfepy.fem.mesh import Mesh
        cmm = Mesh.from_region( control_region, domain.mesh )
        dmm = Mesh.from_region( design_region, domain.mesh )
        cmm.write( 'control.mesh', io = 'auto' )
        dmm.write( 'design.mesh', io = 'auto' )

        SOFC = ShapeOptimFlowCase
        obj = SOFC( sp_boxes = sp_boxes,
                    dsg_vars = dsg_vars,
                    problem_type = opts.problem,
                    objective_function_type = opts.objective_function,
                    var_map = opts.var_map,
                    nu = nu,
                    use_mesh_velocity = opts.use_mesh_velocity,
                    save_dir = opts.save_dir,
                    save_iter_sols = opts.save_iter_sols,
                    save_control_points = opts.save_control_points,
                    save_dsg_vars = opts.save_dsg_vars,
                    test_terms_if_test = opts.test_terms_if_test )

        equations = getattr( conf, '_'.join( ('equations_sensitivity',
                                              opts.problem,
                                              opts.objective_function) ) )

        obj.obj_fun_term = equations['objective']
        obj.sens_terms = equations['sensitivity']

        obj.n_var = dsg_vars.val.shape[0]

        return obj
    from_conf = staticmethod( from_conf )

    ##
    # 07.03.2006, c
    # 12.04.2006
    def generate_mesh_velocity( self, shape, idsgs = None ):
        if self.use_mesh_velocity:
            assert_( shape == self.nu.shape )
            yield self.nu

        else:
            for idsg in idsgs:
#                print 'idsg:', idsg
                yield self.sp_boxes.interp_mesh_velocity( shape, self.dsg_vars,
                                                       idsg )

    ##
    # created: 25.01.2006
    # last revision: 21.12.2007
    def obj_fun( self, vec_dp, conf, domain, variables, materials,
                data = None ):
        """can use both apb, dpb variables"""
        set_state_to_vars( variables, self.var_map, vec_dp )

        val = eva.eval_term( None, self.obj_fun_term, conf, domain,
                             variables, materials, None, **data )
        return nm.squeeze( val )
        
    def sensitivity(self, vec_dp, vec_ap, conf, domain, variables, materials,
                    data=None, select=None):
        """can use both apb, dpb variables"""
        set_state_to_vars( variables, self.var_map, vec_dp, vec_ap )

        dim = self.sp_boxes.dim
        n_mesh_nod = domain.shape.n_nod

        if select is None:
            idsgs = nm.arange( self.dsg_vars.n_dsg, dtype = nm.int32 )
        else:
            idsgs = select

        sa = []

        pbar = MyBar('sensitivity:')
        pbar.init(len(idsgs))

        materials = materials.semideep_copy()
        shape = (n_mesh_nod, dim)
        for ii, nu in enumerate(self.generate_mesh_velocity(shape, idsgs)):
            pbar.update(ii)
            vec_nu = nu.ravel()
            variables['Nu'].data_from_data(vec_nu)
            data.update({'mode' : 1})

            ## from sfepy.base.ioutils import write_vtk
            ## cc = nla.norm( vec_nu )
            ## nun = nu / cc
            ## out = {'v' : Struct( mode = 'vertex', data = nun,
            ##                      ap_name = 'nic', dof_types = (0,1,2) )}
            ## fd = open( 'anim/pert_%03d.pvtk' % (ii+1), 'w' )
            ## write_vtk( fd, domain.mesh, out )
            ## fd.close()
            ## print ii
            
            val = eva.eval_term(None, self.sens_terms, conf, domain,
                                variables, materials, None,
                                copy_materials=False,
                                update_materials=ii==0, **data)

            sa.append( val )

        vec_sa = nm.array( sa, nm.float64 )
        return vec_sa

    ##
    # created:       27.02.2006
    # last revision: 02.01.2008
    def check_custom_sensitivity( self, term_desc, idsg, delta, vec_dp, vec_ap,
                                pb, data = None ):

        variables = pb.variables
        domain = pb.domain
        regions = domain.regions
        materials = pb.materials
        
        # Set admissible state variables.
        set_state_to_vars( variables, self.var_map, vec_dp, vec_ap )

        if data is None:
            data = {}
        else:
            data = copy( data )
            
        dim = self.sp_boxes.dim
        n_mesh_nod = domain.shape.n_nod

        a_grad = []
        d_grad = []

        coors0 = domain.mesh.nod0[:,:-1].copy() 

        for nu in self.generate_mesh_velocity( (n_mesh_nod, dim), [idsg] ):
            vec_nu = nu.ravel()
            variables['Nu'].data_from_data( vec_nu, slice( 0, len( vec_nu ) ) )

            data['mode'] = 1

            aux = eva.eval_term_op( None, term_desc, pb, **data )
            a_grad.append( aux )

            data['mode'] = 0

            coorsp = coors0 + delta * nu
            pb.set_mesh_coors( coorsp, update_state = True )
            valp = eva.eval_term_op( None, term_desc, pb, **data )

            coorsm = coors0 - delta * nu
            pb.set_mesh_coors( coorsm, update_state = True )
            valm = eva.eval_term_op( None, term_desc, pb, **data )
            d_grad.append( 0.5 * (valp - valm) / delta )
            
        pb.set_mesh_coors( coors0, update_state = True )

        a_grad = nm.array( a_grad, nm.float64 )
        d_grad = nm.array( d_grad, nm.float64 )

        output( term_desc + ':' )
        output( '       a: %.8e' % a_grad )
        output( '       d: %.8e' % d_grad )
        output( '-> ratio:', a_grad / d_grad )
        pause()

    ##
    # created:       25.01.2006
    # last revision: 02.01.2008
    def check_sensitivity( self, idsgs, delta, vec_dp, vec_ap, nls,
                          conf, domain, variables, materials, data = None ):
#        print data; pause()
        if data is None:
            data = {}
        else:
            data = copy( data )
        
        dpb = nls.evaluator.problem
        
        a_grad = self.sensitivity( vec_dp, vec_ap, conf, dpb.domain,
                                  dpb.variables, dpb.materials,
                                  data = data, select = idsgs )
        output( a_grad )
        
        d_grad = []

        dim = self.sp_boxes.dim
        n_mesh_nod = domain.shape.n_nod

        coors0 = dpb.get_mesh_coors().copy()
        for nu in self.generate_mesh_velocity( (n_mesh_nod, dim), idsgs ):

            coorsp = coors0 + delta * nu
            dpb.set_mesh_coors( coorsp, update_state = True )

            vec_dp = dpb.create_state_vector()
            dpb.apply_ebc( vec_dp )
            vec_dp = nls( vec_dp )

            valp = self.obj_fun( vec_dp, conf, domain, variables, materials,
                                data )
            output( 'obj_fun+:', valp )

##             out = dm_state_to_output( vec_dp, variables, maps.di, domain )
##             fd = open( trunk + '_sol.vtk', 'w' )
## #            write_vtk( fd, self.char_mesh, out )
##             write_vtk( fd, domain.mesh, out )
##             fd.close()

            coorsm = coors0 - delta * nu
            dpb.set_mesh_coors( coorsm, update_state = True )

            vec_dp = dpb.create_state_vector()
            dpb.apply_ebc( vec_dp )
            vec_dp = nls( vec_dp )

            valm = self.obj_fun( vec_dp, conf, domain, variables, materials,
                                data )
            output( 'obj_fun-:', valm )
            
##             out = dm_state_to_output( vec_dp, variables, maps.di, domain )
##             fd = open( trunk + '_sol.vtk', 'w' )
##             write_vtk( fd, domain.mesh, out )
##             fd.close()
## #            pause()

            d_grad.append( 0.5 * (valp - valm) / delta )
            
        ##
        # Restore original mesh coordinates.
        dpb.set_mesh_coors( coors0, update_state = True )

        d_grad = nm.array( d_grad, nm.float64 )

        output( 'a: %.8e' % a_grad )
        output( 'd: %.8e' % d_grad )
        output( a_grad / d_grad )
        pause()
