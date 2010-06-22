import os.path as op

from sfepy.base.base import *
from sfepy.base.progressbar import MyBar
import sfepy.fem.evaluate as eva
import freeFormDef as ffd
from sfepy.terms import CharacteristicFunction

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
                             var_data=None,
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
            
        # Restore materials.
        pb.time_update()

        # Solve direct/adjoint problem.
        vec = problem.solve(var_data=var_data)

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

    val = shape_opt.obj_fun(vec_dp, dpb, data=data)
#    print 'OF ok'
    return val

##
# created:       19.04.2006
# last revision: 21.12.2007
def obj_fun_grad( design, shape_opt, dpb, apb, opts ):

#    print 'OFG'
    data = {}
    vec_dp = solve_problem_for_design(dpb, design, shape_opt, opts )

    var_data = dpb.equations.get_state_parts(vec_dp)
    var_data = remap_dict(var_data, shape_opt.var_map)

    vec_ap = solve_problem_for_design(apb, design, shape_opt, opts,
                                      var_data=var_data,
                                      use_cache=False, is_mesh_update=False)

    vec_sa = shape_opt.sensitivity(var_data, vec_ap, apb, data=data)
#    print 'OFG ok'
    return vec_sa

##
# 13.04.2006, c
# 26.07.2006
# 20.03.2007
# 22.10.2007
# 24.10.2007
# 25.10.2007
def test_terms(idsgs, delta, shape_opt, dp_var_data, vec_ap, pb):

    dd = {}
    ccs = shape_opt.check_custom_sensitivity

    ccs('d_sd_st_grad_div.i2.Omega_D( stabil.gamma, w, u, Nu, mode )',
        idsgs, delta, dp_var_data, vec_ap, pb, dd)

    ccs('d_sd_st_supg_c.i1.Omega_D( stabil.delta, w, w, u, Nu, mode )',
        idsgs, delta, dp_var_data, vec_ap, pb, dd)

    ccs('d_sd_st_pspg_c.i1.Omega_D( stabil.tau, r, w, u, Nu, mode )',
        idsgs, delta, dp_var_data, vec_ap, pb, dd)

    ccs('d_sd_st_pspg_p.i1.Omega( stabil.tau, r, p, Nu, mode )',
        idsgs, delta, dp_var_data, vec_ap, pb, dd)
    ccs('d_sd_st_pspg_p.i1.Omega( stabil.tau, r, r, Nu, mode )',
        idsgs, delta, dp_var_data, vec_ap, pb, dd)
    ccs('d_sd_st_pspg_p.i1.Omega( stabil.tau, p, p, Nu, mode )',
        idsgs, delta, dp_var_data, vec_ap, pb, dd)

    ccs('d_sd_test_pq.i1.Omega_D( p, r, Nu, mode )',
        idsgs, delta, dp_var_data, vec_ap, pb, dd)

    ccs('d_sd_div.i1.Omega_D( u, r, Nu, mode )',
        idsgs, delta, dp_var_data, vec_ap, pb, dd)

    ccs('d_sd_div.i1.Omega_D( w, p, Nu, mode )',
        idsgs, delta, dp_var_data, vec_ap, pb, dd)

    ccs('d_sd_div_grad.i2.Omega_D( one.val, fluid.viscosity, u, w, Nu, mode )',
        idsgs, delta, dp_var_data, vec_ap, pb, dd)

    ccs('d_sd_convect.i2.Omega_D( u, w, Nu, mode )',
        idsgs, delta, dp_var_data, vec_ap, pb, dd)

    # Restore materials.
    pb.time_update()

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

    def obj_fun(self, vec_dp, dpb, data=None):
        """
        Objective function evaluation for given direct problem state.
        """
        var_data = dpb.equations.get_state_parts(vec_dp)
        var_data = remap_dict(var_data, self.var_map)

        val = dpb.evaluate(self.obj_fun_term, copy_materials=False,
                           var_names=self.var_map.keys(), **var_data)

        return nm.squeeze( val )
        
    def sensitivity(self, dp_var_data, vec_ap, apb, data=None, select=None):
        """
        Sensitivity of objective function evaluation for given direct
        and adjoint problem states.
        """
        var_data = apb.equations.get_state_parts(vec_ap)
        var_data.update(dp_var_data)
        var_names = var_data.keys() + ['Nu']

        var_data.update(data)
        
        dim = self.sp_boxes.dim
        n_mesh_nod = apb.domain.shape.n_nod

        if select is None:
            idsgs = nm.arange( self.dsg_vars.n_dsg, dtype = nm.int32 )
        else:
            idsgs = select

        sa = []

        pbar = MyBar('sensitivity:')
        pbar.init(len(idsgs))

        shape = (n_mesh_nod, dim)
        for ii, nu in enumerate(self.generate_mesh_velocity(shape, idsgs)):
            pbar.update(ii)
            var_data['Nu'] = nu.ravel()
            var_data['mode'] = 1

            ## from sfepy.base.ioutils import write_vtk
            ## cc = nla.norm( vec_nu )
            ## nun = nu / cc
            ## out = {'v' : Struct( mode = 'vertex', data = nun,
            ##                      ap_name = 'nic', dof_types = (0,1,2) )}
            ## fd = open( 'anim/pert_%03d.pvtk' % (ii+1), 'w' )
            ## write_vtk( fd, domain.mesh, out )
            ## fd.close()
            ## print ii

            val = apb.evaluate(self.sens_terms, copy_materials=False,
                               update_materials=ii==0,
                               new_geometries=ii==0,
                               var_names=var_names, **var_data)

            sa.append( val )

        vec_sa = nm.array( sa, nm.float64 )
        return vec_sa

    def check_custom_sensitivity(self, term_desc, idsg, delta,
                                 dp_var_data, vec_ap, pb, data=None):

        domain = pb.domain
        regions = domain.regions
        materials = pb.materials

        if data is None:
            data = {}
        else:
            data = copy( data )

        var_data = pb.equations.get_state_parts(vec_ap)
        var_data.update(dp_var_data)
        var_names = var_data.keys() + ['Nu']

        var_data.update(data)

        dim = self.sp_boxes.dim
        n_mesh_nod = domain.shape.n_nod

        a_grad = []
        d_grad = []

        coors0 = domain.mesh.coors

        for nu in self.generate_mesh_velocity( (n_mesh_nod, dim), [idsg] ):
            var_data['Nu'] = nu.ravel()
            var_data['mode'] = 1

            aux = pb.evaluate(term_desc, copy_materials=False,
                              update_materials=True,
                              var_names=var_names, **var_data)
            a_grad.append( aux )

            var_data['mode'] = 0

            coorsp = coors0 + delta * nu
            pb.set_mesh_coors( coorsp, update_state = True )
            valp = pb.evaluate(term_desc, copy_materials=False,
                               update_materials=False,
                               var_names=var_names, **var_data)

            coorsm = coors0 - delta * nu
            pb.set_mesh_coors( coorsm, update_state = True )
            valm = pb.evaluate(term_desc, copy_materials=False,
                               update_materials=False,
                               var_names=var_names, **var_data)

            d_grad.append( 0.5 * (valp - valm) / delta )
            
        pb.set_mesh_coors( coors0, update_state = True )

        a_grad = nm.array( a_grad, nm.float64 )
        d_grad = nm.array( d_grad, nm.float64 )

        output( term_desc + ':' )
        output( '       a: %.8e' % a_grad )
        output( '       d: %.8e' % d_grad )
        output( '-> ratio:', a_grad / d_grad )
        pause()

    def check_sensitivity(self, idsgs, delta, dp_var_data, vec_ap, dpb,
                          apb, data=None):
        if data is None:
            data = {}
        else:
            data = copy( data )
        
        a_grad = self.sensitivity(dp_var_data, vec_ap, apb,
                                  data=data, select=idsgs)
        output( a_grad )
        
        d_grad = []

        dim = self.sp_boxes.dim
        n_mesh_nod = dpb.domain.shape.n_nod

        coors0 = dpb.get_mesh_coors().copy()
        for nu in self.generate_mesh_velocity( (n_mesh_nod, dim), idsgs ):

            coorsp = coors0 + delta * nu
            dpb.set_mesh_coors( coorsp, update_state = True )

            dpb.time_update()
            vec_dp = dpb.solve()

            valp = self.obj_fun(vec_dp, dpb, data)
            output( 'obj_fun+:', valp )

            coorsm = coors0 - delta * nu
            dpb.set_mesh_coors( coorsm, update_state = True )

            dpb.time_update()
            vec_dp = dpb.solve()

            valm = self.obj_fun(vec_dp, dpb, data)
            output( 'obj_fun-:', valm )
            
            d_grad.append( 0.5 * (valp - valm) / delta )
            
        ##
        # Restore original mesh coordinates.
        dpb.set_mesh_coors( coors0, update_state = True )

        d_grad = nm.array( d_grad, nm.float64 )

        output( 'a: %.8e' % a_grad )
        output( 'd: %.8e' % d_grad )
        output( a_grad / d_grad )
        pause()
