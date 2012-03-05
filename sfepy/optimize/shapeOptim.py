import os.path as op

import numpy as nm

from sfepy.base.base import output, assert_, remap_dict, pause, Struct
from sfepy.base.progressbar import MyBar
from sfepy.fem.equations import get_expression_arg_names
from sfepy.fem.evaluate import eval_equations
import freeFormDef as ffd

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
    pb.set_mesh_coors( coors, update_fields=False )

    pb.domain.mesh.write( op.join( shape_opt.save_dir, 'design.%03d.mesh' )\
                          % shape_opt.cache.i_mesh, io = 'auto' )
    shape_opt.cache.i_mesh += 1
    try:
        pb.set_mesh_coors( coors, update_fields=True )
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
        state = shape_opt.cache.state
    else:
        if is_mesh_update:
            ok = update_mesh( shape_opt, pb, design )
            if not ok: return None

        # Solve direct/adjoint problem.
        state = problem.solve(var_data=var_data)

        if use_cache:
            shape_opt.cache.state = state.copy()
            shape_opt.cache.design = design.copy()

            if shape_opt.save_iter_sols:
                pb.save_state( op.join( shape_opt.save_dir, 'direct.%03d.vtk' )\
                              % (shape_opt.cache.i_mesh - 1), state )

            if shape_opt.save_control_points:
                aux = shape_opt.sp_boxes.create_mesh_from_control_points()
                aux.write( op.join( shape_opt.save_dir, 'cp.%03d.vtk' )\
                           % (shape_opt.cache.i_mesh - 1), io = 'auto' )

            if shape_opt.save_dsg_vars:
                filename = op.join( shape_opt.save_dir, 'dsg.%03d.txt' )\
                           % (shape_opt.cache.i_mesh - 1)
                shape_opt.dsg_vars.val.tofile( filename, ' ' )

    return state

def obj_fun(design, shape_opt, opts):
    """
    The objective function evaluation.
    """
    state_dp = solve_problem_for_design(shape_opt.dpb, design, shape_opt, opts)
    if state_dp is None:
        return None

    val = shape_opt.obj_fun(state_dp)

    return val

def obj_fun_grad(design, shape_opt, opts):
    """
    The objective function gradient evaluation.
    """
    state_dp = solve_problem_for_design(shape_opt.dpb, design, shape_opt, opts)

    var_data = state_dp.get_parts()
    var_data = remap_dict(var_data, shape_opt.var_map)

    state_ap = solve_problem_for_design(shape_opt.apb, design, shape_opt, opts,
                                        var_data=var_data,
                                        use_cache=False, is_mesh_update=False)

    vec_sa = shape_opt.sensitivity(var_data, state_ap)

    return vec_sa

def test_terms(idsgs, delta, shape_opt, dp_var_data, state_ap):
    """
    Test individual shape derivative terms.
    """

    def ccs(*args):
        try:
            shape_opt.check_custom_sensitivity(*args)

        except:
            output('running test failed!')

    ccs('d_sd_st_grad_div.i2.Omega_D( stabil.gamma, w, u, Nu )',
        idsgs, delta, dp_var_data, state_ap)

    ccs('d_sd_st_supg_c.i1.Omega_D( stabil.delta, w, w, u, Nu )',
        idsgs, delta, dp_var_data, state_ap)

    ccs('d_sd_st_pspg_c.i1.Omega_D( stabil.tau, r, w, u, Nu )',
        idsgs, delta, dp_var_data, state_ap)

    ccs('d_sd_st_pspg_p.i1.Omega_D( stabil.tau, r, p, Nu )',
        idsgs, delta, dp_var_data, state_ap)
    ccs('d_sd_st_pspg_p.i1.Omega_D( stabil.tau, r, r, Nu )',
        idsgs, delta, dp_var_data, state_ap)
    ccs('d_sd_st_pspg_p.i1.Omega_D( stabil.tau, p, p, Nu )',
        idsgs, delta, dp_var_data, state_ap)

    ccs('d_sd_dot_scalar.i1.Omega_D( p, r, Nu )',
        idsgs, delta, dp_var_data, state_ap)

    ccs('d_sd_div.i1.Omega_D( u, r, Nu )',
        idsgs, delta, dp_var_data, state_ap)
    ccs('d_sd_div.i1.Omega_D( w, p, Nu )',
        idsgs, delta, dp_var_data, state_ap)

    ccs('d_sd_div_grad.i2.Omega_D( one.val, fluid.viscosity, u, w, Nu )',
        idsgs, delta, dp_var_data, state_ap)

    ccs('d_sd_convect.i2.Omega_D( u, w, Nu )',
        idsgs, delta, dp_var_data, state_ap)

##
# 25.01.2006, c
class ShapeOptimFlowCase( Struct ):

    def from_conf(conf, dpb, apb):

        opts = conf.options
        regions = dpb.domain.regions
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
        cmm = Mesh.from_region(control_region, dpb.domain.mesh)
        dmm = Mesh.from_region(design_region, dpb.domain.mesh)
        cmm.write( 'control.mesh', io = 'auto' )
        dmm.write( 'design.mesh', io = 'auto' )

        SOFC = ShapeOptimFlowCase
        obj = SOFC(dpb=dpb, apb=apb,
                   sp_boxes=sp_boxes,
                   dsg_vars=dsg_vars,
                   problem_type=opts.problem,
                   objective_function_type=opts.objective_function,
                   var_map=opts.var_map,
                   nu=nu,
                   use_mesh_velocity=opts.use_mesh_velocity,
                   save_dir=opts.save_dir,
                    save_iter_sols=opts.save_iter_sols,
                   save_control_points=opts.save_control_points,
                   save_dsg_vars=opts.save_dsg_vars,
                   test_terms_if_test=opts.test_terms_if_test)

        equations = getattr( conf, '_'.join( ('equations_sensitivity',
                                              opts.problem,
                                              opts.objective_function) ) )

        obj.obj_fun_term = equations['objective']
        obj.sens_terms = equations['sensitivity']

        obj.n_var = dsg_vars.val.shape[0]

        obj.create_evaluables()

        return obj
    from_conf = staticmethod( from_conf )

    def create_evaluables(self):
        variables = self.dpb.create_variables().as_dict()

        possible_mat_names = get_expression_arg_names(self.obj_fun_term)
        materials = self.dpb.create_materials(possible_mat_names).as_dict()

        aux = self.dpb.create_evaluable(self.obj_fun_term,
                                        try_equations=False,
                                        var_dict=variables,
                                        **materials)
        self.of_equations, self.of_variables = aux

        possible_mat_names = get_expression_arg_names(self.sens_terms)
        materials = self.dpb.create_materials(possible_mat_names).as_dict()

        aux = self.apb.create_evaluable(self.sens_terms,
                                        try_equations=False,
                                        var_dict=variables,
                                        **materials)
        self.ofg_equations, self.ofg_variables = aux

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

    def obj_fun(self, state_dp):
        """
        Objective function evaluation for given direct problem state.
        """
        var_data = state_dp.get_parts()
        var_data = remap_dict(var_data, self.var_map)

        self.of_equations.set_data(var_data, ignore_unknown=True)

        val = eval_equations(self.of_equations, self.of_variables)

        return nm.squeeze( val )

    def sensitivity(self, dp_var_data, state_ap, select=None):
        """
        Sensitivity of objective function evaluation for given direct
        and adjoint problem states.
        """
        apb = self.apb

        var_data = state_ap.get_parts()
        var_data.update(dp_var_data)

        self.ofg_equations.set_data(var_data, ignore_unknown=True)

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
            self.ofg_variables['Nu'].data_from_any(nu.ravel())

            ## from sfepy.base.ioutils import write_vtk
            ## cc = nla.norm( vec_nu )
            ## nun = nu / cc
            ## out = {'v' : Struct( mode = 'vertex', data = nun,
            ##                      ap_name = 'nic', dof_types = (0,1,2) )}
            ## fd = open( 'anim/pert_%03d.pvtk' % (ii+1), 'w' )
            ## write_vtk( fd, domain.mesh, out )
            ## fd.close()
            ## print ii

            val = eval_equations(self.ofg_equations, self.ofg_variables,
                                 term_mode=1)

            sa.append( val )

        vec_sa = nm.array( sa, nm.float64 )
        return vec_sa

    def check_custom_sensitivity(self, term_desc, idsg, delta,
                                 dp_var_data, state_ap):
        pb = self.apb

        domain = pb.domain

        possible_mat_names = get_expression_arg_names(term_desc)
        materials = self.dpb.create_materials(possible_mat_names).as_dict()

        variables = self.ofg_equations.variables
        aux = self.dpb.create_evaluable(term_desc,
                                        try_equations=False,
                                        var_dict=variables,
                                        **materials)
        check0_equations, check0_variables = aux

        aux = self.dpb.create_evaluable(term_desc,
                                        try_equations=False,
                                        var_dict=variables,
                                        **materials)
        check1_equations, check1_variables = aux

        var_data = state_ap.get_parts()
        var_data.update(dp_var_data)

        check0_equations.set_data(var_data, ignore_unknown=True)
        check1_equations.set_data(var_data, ignore_unknown=True)

        dim = self.sp_boxes.dim
        n_mesh_nod = domain.shape.n_nod

        a_grad = []
        d_grad = []

        coors0 = domain.mesh.coors

        for nu in self.generate_mesh_velocity( (n_mesh_nod, dim), [idsg] ):
            check1_variables['Nu'].data_from_any(nu.ravel())

            aux = eval_equations(check1_equations, check1_variables,
                                 term_mode=1)
            a_grad.append( aux )

            coorsp = coors0 + delta * nu
            pb.set_mesh_coors( coorsp, update_fields=True )
            valp = eval_equations(check0_equations, check0_variables,
                                  term_mode=0)

            coorsm = coors0 - delta * nu
            pb.set_mesh_coors( coorsm, update_fields=True )
            valm = eval_equations(check0_equations, check0_variables,
                                  term_mode=0)

            d_grad.append( 0.5 * (valp - valm) / delta )

        pb.set_mesh_coors( coors0, update_fields=True )

        a_grad = nm.array( a_grad, nm.float64 )
        d_grad = nm.array( d_grad, nm.float64 )

        output( term_desc + ':' )
        output( '       a: %.8e' % a_grad )
        output( '       d: %.8e' % d_grad )
        output( '-> ratio:', a_grad / d_grad )
        pause()

    def check_sensitivity(self, idsgs, delta, dp_var_data, state_ap):
        a_grad = self.sensitivity(dp_var_data, state_ap, select=idsgs)
        output( a_grad )

        d_grad = []

        dpb = self.dpb
        dim = self.sp_boxes.dim
        n_mesh_nod = dpb.domain.shape.n_nod

        coors0 = dpb.get_mesh_coors().copy()
        for nu in self.generate_mesh_velocity( (n_mesh_nod, dim), idsgs ):

            coorsp = coors0 + delta * nu
            dpb.set_mesh_coors( coorsp, update_fields=True )

            dpb.time_update()
            state_dp = dpb.solve()

            valp = self.obj_fun(state_dp)
            output( 'obj_fun+:', valp )

            coorsm = coors0 - delta * nu
            dpb.set_mesh_coors( coorsm, update_fields=True )

            dpb.time_update()
            state_dp = dpb.solve()

            valm = self.obj_fun(state_dp)
            output( 'obj_fun-:', valm )

            d_grad.append( 0.5 * (valp - valm) / delta )

        ##
        # Restore original mesh coordinates.
        dpb.set_mesh_coors( coors0, update_fields=True )

        d_grad = nm.array( d_grad, nm.float64 )

        output( 'a: %.8e' % a_grad )
        output( 'd: %.8e' % d_grad )
        output( a_grad / d_grad )
        pause()
