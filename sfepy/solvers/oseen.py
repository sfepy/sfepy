from sfepy.base.base import *
from sfepy.solvers.solvers import NonlinearSolver
from nls import conv_test
import sfepy.base.plotutils as plu

##
# 26.07.2007, c
def are_close( a, b, rtol = 0.2, atol = 1e-8 ):
    return False
#    return abs( a - b ) <= max( atol, rtol * abs( b ) )

##
# 26.07.2007, c
def scale_matrix( mtx, indx, factor ):
    ptr0 = mtx.indptr[indx.start]
    ptr1 = mtx.indptr[indx.stop]
    mtx.data[ptr0:ptr1] *= factor

_dimater_modes = {'edge' : 0, 'volume' : 1, 'max' : 2}

##
# c: 01.08.2007, r: 15.01.2008
def create_stabil_data( problem, fluid_name, stabil_name, eq_name1, eq_name2 ):

    ns = {}
    term = problem.equations[eq_name1].terms['dw_lin_convect']

    ns['fluid'] = 'fluid'
    ns['v'] = term.get_virtual_name()
    ns['b'] = term.get_parameter_names()[0]
    ns['u'] = term.get_state_names()[0]
    ns['omega'] = term.region.name
    
    term = problem.equations[eq_name1].terms['dw_stokes']
    ns['p'] = term.get_state_names()[0]

    term = problem.equations[eq_name2].terms['dw_stokes']
    ns['q'] = term.get_virtual_name()

    ii = {}
    ii['u'] = problem.variables.get_indx( ns['u'] )
    ii['us'] = problem.variables.get_indx( ns['u'], stripped = True )
    ii['ps'] = problem.variables.get_indx( ns['p'], stripped = True )

    stabil = problem.materials[stabil_name]
    mat = problem.materials[ns['fluid']]

    viscosity = mat.viscosity

    c_friedrichs = problem.domain.get_diameter()
    sigma = 1e-12 # 1 / dt.

#    print c_friedrichs

    def mat_fun( ts, coor, region, ig, b_norm = 1.0, fixed_data = None ):
        if fixed_data is not None:
            return fixed_data[ig]

        print '|b|_max (mat_fun):', b_norm
        gamma = viscosity + b_norm * c_friedrichs

        data = {}
        if stabil.gamma is None:
            data['gamma'] = stabil.gamma_mul * gamma
        else:
            data['gamma'] = nm.asarray( stabil.gamma_mul * stabil.gamma,
                                        dtype = nm.float64 )
        

        if stabil.delta is None:
            term = problem.equations[eq_name1].terms['dw_lin_convect']
            for ig in term.iter_groups():
                # This sets term.ig - for 1 group only!!!
                break
            var = problem.variables[ns['u']]
            ap, vg = var.get_approximation( term.get_current_group(), 'Volume' )
            delta = 1.0
            mode = _dimater_modes[stabil.diameter_mode]
            cells = stabil.region.get_cells( ig )
            diameters2 = problem.domain.get_element_diameters( ig, cells, vg,
                                                             mode )
            val1 = min( 1.0, 1.0 / sigma )
            val2 = sigma * c_friedrichs**2
            val3 = (b_norm**2) * min( (c_friedrichs**2) / viscosity, 1.0 / sigma )
#            print val1, gamma, val2, val3
            delta = stabil.delta_mul * val1 * diameters2 / (gamma + val2 + val3)
            data['diameters2'] = diameters2
            data['delta'] = delta
        else:
            val = stabil.delta_mul * stabil.delta
            for ii in range( len( stabil.igs ) ):
                data['delta'] = nm.asarray( val, dtype = nm.float64 )
        
        if stabil.tau is None:
            data['tau'] = stabil.tau_red * data['delta']
        else:
            data['tau'] = nm.asarray( stabil.tau_mul * stabil.tau,
                                      dtype = nm.float64 )

        return data

    stabil.set_function( mat_fun )

    return stabil, ns, ii

##
# 11.10.2007, c
class Oseen( NonlinearSolver ):
    name = 'nls.oseen'

    def process_conf( conf ):
        """
        Missing items are set to default values.
        
        Example configuration, all items:
        
        solver_1 = {
            'name' : 'oseen',
            'kind' : 'nls.oseen',

            'adimensionalize' : False,
            'check_navier_stokes_rezidual' : False,

            'fluid_mat_name' : 'fluid',
            'stabil_mat_name' : 'stabil',
            'lin_convect_eq_name' : 'balance',
            'div_eq_name' : 'incompressibility',

            'i_max'      : 10,
            'eps_a'      : 1e-8,
            'eps_r'      : 1.0,
            'macheps'   : 1e-16,
            'lin_red'    : 1e-2, # Linear system error < (eps_a * lin_red).
            'is_plot'    : False,
        }
        """
        get = conf.get_default_attr

        # Compulsory.
        fluid_mat_name = get( 'fluid_mat_name', None,
                              'missing "fluid_mat_name" in options!' )
        stabil_mat_name = get( 'stabil_mat_name', None,
                               'missing "stabil_mat_name" in options!' )
        lin_convect_eq_name = get( 'lin_convect_eq_name', None,
                                   'missing "lin_convect_eq_name" in options!' )
        div_eq_name = get( 'div_eq_name', None,
                           'missing "div_eq_name" in options!' )

        # With defaults.
        adimensionalize = get( 'adimensionalize', False )
        check_navier_stokes_rezidual = get( 'check_navier_stokes_rezidual',
                                            False )
        i_max = get( 'i_max', 1 )
        eps_a = get( 'eps_a', 1e-10 )
        eps_r = get( 'eps_r', 1.0 )
        macheps = get( 'macheps', nm.finfo( nm.float64 ).eps )
        lin_red = get( 'lin_red', 1.0 )
        is_plot = get( 'is_plot', False )

        common = NonlinearSolver.process_conf( conf )
        return Struct( **locals() ) + common
    process_conf = staticmethod( process_conf )

    ##
    # 10.10.2007, c
    def __init__( self, conf, **kwargs ):
        NonlinearSolver.__init__( self, conf, **kwargs )

    ##
    # 26.07.2007, c
    # 31.07.2007
    # 01.08.2007
    # 11.10.2007, from oseen()
    # 29.10.2007
    # 30.10.2007
    # 31.10.2007
    def __call__( self, vec_x0, conf = None, evaluator = None,
                  lin_solver = None, status = None ):
        """"""
        conf = get_default( conf, self.conf )
        evaluator = get_default( evaluator, self.evaluator )
        lin_solver = get_default( lin_solver, self.lin_solver )
        status = get_default( status, self.status )

        if hasattr( conf, 'fixed_data' ):
            fixed_data = conf.fixed_data
        else:
            fixed_data = None

        time_stats = {}

        problem = evaluator.problem

        stabil, ns, ii = create_stabil_data( problem, conf.fluid_mat_name,
                                           conf.stabil_mat_name,
                                           conf.lin_convect_eq_name,
                                           conf.div_eq_name )
        update_var = problem.variables.non_state_data_from_state

        print 'problem size:'
        print '    velocity: %s' % ii['us']
        print '    pressure: %s' % ii['ps']

        vec_x = vec_x0.copy()
        vec_x_prev = vec_x0.copy()
        vec_dx = None

        err0 = -1.0
        it = 0
        while 1:
            update_var( ns['b'], vec_x_prev, ns['u'] )
            vec_b = vec_x_prev[ii['u']]
            b_norm = nla.norm( vec_b, nm.inf )
            print '|b|_max: %.12e' % b_norm

            vec_u = vec_x[ii['u']]
            u_norm = nla.norm( vec_u, nm.inf )
            print '|u|_max: %.2e' % u_norm

            stabil.time_update( None, None, problem.domain,
                               b_norm = b_norm, fixed_data = fixed_data )
            max_pars = stabil.reduce_on_datas( lambda a, b: max( a, b.max() ) )
            print 'stabilization parameters:'
            print '                   gamma: %.12e' % max_pars['gamma']
            print '            max( delta ): %.12e' % max_pars['delta']
            print '              max( tau ): %.12e' % max_pars['tau']
            try:
                print '              max( h^2 ): %.12e' % max_pars['diameters2']
            except:
                pass

            if (not are_close( b_norm, 1.0 )) and conf.adimensionalize:
                adimensionalize = True
            else:
                adimensionalize = False

            tt = time.clock()
            vec_r, ret = evaluator.eval_residual( vec_x )
            time_stats['rezidual'] = time.clock() - tt
            if ret == 0: # OK.
                err = nla.norm( vec_r )
                if it == 0:
                    err0 = err;
                else:
                    err += nla.norm( vec_dx )
            else: # Failure.
                print 'rezidual computation failed for iter %d!' % it
                raise RuntimeError, 'giving up...'

            condition = conv_test( conf, it, err, err0 )
            if condition >= 0:
                break

            if adimensionalize:
                print 'adimensionalizing'
                mat.viscosity = viscosity / b_norm
                vec_r[indx_us] /= b_norm

            tt = time.clock()
            mtx_a, ret = evaluator.eval_tangent_matrix( vec_x )
            time_stats['matrix'] = time.clock() - tt
            if ret != 0:
                raise RuntimeError, 'giving up...'

            tt = time.clock() 
            vec_dx = lin_solver( vec_r, mtx = mtx_a )
            time_stats['solve'] = time.clock() - tt

            vec_e = mtx_a * vec_dx - vec_r
            lerr = nla.norm( vec_e )
            if lerr > (conf.eps_a * conf.lin_red):
                print 'linear system not solved! (err = %e)' % lerr
    #            raise RuntimeError, 'linear system not solved! (err = %e)' % lerr

            if adimensionalize:
                print 'restoring pressure...'
                vec_dx[indx_ps] *= b_norm

            dx_norm = nla.norm( vec_dx )
            print '||dx||: %.2e' % dx_norm

            for kv in time_stats.iteritems():
                print '%10s: %7.2f [s]' % kv


            vec_x_prev = vec_x.copy()
            evaluator.update_vec( vec_x, vec_dx )

            if conf.is_plot:
                plu.pylab.ion()
                plu.pylab.gcf().clear()
                plu.pylab.subplot( 2, 2, 1 )
                plu.pylab.plot( vec_x_prev )
                plu.pylab.ylabel( r'$x_{i-1}$' )
                plu.pylab.subplot( 2, 2, 2 )
                plu.pylab.plot( vec_r )
                plu.pylab.ylabel( r'$r$' )
                plu.pylab.subplot( 2, 2, 4 )
                plu.pylab.plot( vec_dx )
                plu.pylab.ylabel( r'$\_delta x$' )
                plu.pylab.subplot( 2, 2, 3 )
                plu.pylab.plot( vec_x )
                plu.pylab.ylabel( r'$x_i$' )
                plu.pylab.draw()
                plu.pylab.ioff()
                pause()

            it += 1

        if conf.check_navier_stokes_rezidual:
    ##         update_var( b_name, vec_x_prev, u_name )
    ## #        update_var( b_name, vec_x, u_name )
    ##         vec_rns1, ret = residual( vec_x, context )
    ##         err_ns = nla.norm( vec_rns1 )
    ##         print '"Oseen" rezidual: %.8e' % err_ns


            t1 = '+ dw_div_grad.%s( %s, %s, %s )' % (ns['omega'],
                                                     ns['fluid'],
                                                     ns['v'], ns['u'])
##             t2 = '+ dw_lin_convect.%s( %s, %s, %s )' % (ns['omega'],
##                                                         ns['v'], b_name, ns['u'])
            t2 = '+ dw_convect.%s( %s, %s )' % (ns['omega'], ns['v'], ns['u'])
            t3 = '- dw_grad.%s( %s, %s )' % (ns['omega'], ns['v'], ns['p'])
            t4 = 'dw_div.%s( %s, %s )' % (ns['omega'], ns['q'], ns['u'])
            equations = {
                'balance' : ' '.join( (t1, t2, t3) ),
                'incompressibility' : t4,
            }
            problem.set_equations( equations )
            vec_rns0, ret = evaluator.eval_residual( vec_x0 )
            vec_rns, ret = evaluator.eval_residual( vec_x )
            if ret:
                print 'Navier-Stokes rezidual computation failed!'
                err_ns = err_ns0 = None
            else:
                err_ns0 = nla.norm( vec_rns0 )
                err_ns = nla.norm( vec_rns )
            print 'Navier-Stokes rezidual0: %.8e' % err_ns0
            print 'Navier-Stokes rezidual : %.8e' % err_ns
            print 'b - u: %.8e' % nla.norm( vec_b - vec_u )
            print condition
    ##         print vec_rns - vec_rns1
            plu.pylab.ion()
            plu.pylab.gcf().clear()
            plu.pylab.plot( vec_rns )
    ##         plu.pylab.gcf().clear()
    ##         plu.pylab.plot( vec_rns1 )
            plu.pylab.draw()
            plu.pylab.ioff()
            pause()
        else:
            err_ns = None

        if status is not None:
            status['time_stats'] = time_stats
            status['err0'] = err0
            status['err'] = err
            status['err_ns'] = err_ns
            status['condition'] = condition

        return vec_x
