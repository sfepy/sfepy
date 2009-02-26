#!/usr/bin/env python
# 11.07.2006, c 
import os.path as op
import shutil
from optparse import OptionParser

import init_sfepy
from sfepy.base.base import *
from sfepy.base.conf import ProblemConf
from sfepy.fem.problemDef import ProblemDefinition
from sfepy.fem.evaluate import eval_term_op
import sfepy.homogenization.pfdpm as pfdpm
from sfepy.homogenization.coefs_base import MiniAppBase
from sfepy.homogenization.engine import HomogenizationEngine
from sfepy.applications import SimpleApp
from sfepy.base.conf import get_standard_keywords
from sfepy.solvers.generic import solve_evolutionary_op, solve_stationary_op

class Volume( MiniAppBase ):

    def __call__( self, problem = None ):
        problem = get_default( problem, self.problem )
        problem.select_variables( self.variables )

        volume = eval_term_op( None, self.expression, problem )

        return volume

def get_evp( key, cache_evp, problem, conf, equivalence = None ):
    """
    Calls pfdpm.solve_pressure_eigenproblem() if the same or equivalent EVP
    was not already solved - in that case returns the cached EVP.
    """
    evp = None
    if equivalence is None:
        if cache_evp.has_key( key ):
            evp = cache_evp[key]
    else:
        for key2 in equivalence[key]:
            if key2 in cache_evp:
                evp = cache_evp[key2]
                cache_evp[key] = evp
                break
                
    if evp is None:
        ebc = pfdpm.select_by_names( conf.ebcs, conf.ebc_sets[key] )
        epbc = pfdpm.select_by_names( conf.epbcs, conf.epbc_sets[key] )

        solve = pfdpm.get_matrix_parts
        matrices = solve( problem, ebc, epbc, conf.variables,
                          conf.equations_time )
        solve = pfdpm.solve_pressure_eigenproblem
        evp = solve( matrices, conf.options.eig_problem,
                     conf.options.n_eigs,
                     conf.options.check.get( 'diagonalization', False ) )
        evp.ebc, evp.epbc = ebc, epbc
        cache_evp[key] = evp
    else:
        # Just create/restore equation mappings.
        variables = select_by_names( conf.variables, ['uc', 'vc', 'pc', 'qc'] )
        problem.set_variables( variables )
        problem.set_equations( conf.equations_time )
        problem.time_update( conf_ebc = evp.ebc, conf_epbc = evp.epbc )

    return evp

def build_evp_equivalence( equivs_in ):
    if equivs_in is None:
        return None

    out = {}
    for equiv in equivs_in:
        for item in equiv:
            if item in out:
                out[item].add( equiv )
            else:
                out[item] = set( equiv )
    return out

class PorousMediaHomogenizationApp( SimpleApp ):

    def process_options( options ):
        """Application options setup. Sets default values for missing
        non-compulsory options."""
        get = options.get_default_attr
        
        print_digits = get( 'print_digits', 3 )

        float_format = get( 'float_format', '%8.3e' )
        coef_save_name = get( 'coef_save_name', 'coefs' )

        coefs = get( 'coefs', None, 'missing "coefs" in options!' )
        requirements = get( 'requirements', None,
                            'missing "requirements" in options!' )
        volume = get( 'volume', None, 'missing "volume" in options!' )

        return Struct( **locals() )
    process_options = staticmethod( process_options )

    def __init__( self, conf, options, output_prefix, **kwargs ):
        SimpleApp.__init__( self, conf, options, output_prefix,
                            init_equations = False )

        self.setup_options()
        self.cached_coefs = None

        output_dir = self.problem.output_dir
        shutil.copyfile( conf._filename,
                         op.join( output_dir, op.basename( conf._filename ) ) )

    def setup_options( self ):
        SimpleApp.setup_options( self )
        po = PorousMediaHomogenizationApp.process_options
        self.app_options += po( self.conf.options )
    
    def call( self, ret_all = False ):
        opts = self.app_options
        
        volume = Volume( 'volume', self.problem, opts.volume )()
        output( 'volume: %.2f' % volume )
        
        he = HomogenizationEngine( self.problem, self.options, volume = volume )

        aux = he( ret_all = ret_all)
        if ret_all:
            coefs, dependencies = aux
        else:
            coefs = aux

        coefs = pfdpm.Coefficients( **coefs.to_dict() )
        coefs.volume = volume
        
        prec = nm.get_printoptions()[ 'precision']
        if hasattr( opts, 'print_digits' ):
            nm.set_printoptions( precision = opts.print_digits )
        print coefs
        nm.set_printoptions( precision = prec )
##        pause()

        coef_save_name = op.join( opts.output_dir, opts.coef_save_name )
        coefs.to_file_hdf5( coef_save_name + '.h5' )
        coefs.to_file_txt( coef_save_name + '.txt',
                           self.conf.options.tex_names,
                           opts.float_format )

        if ret_all:
            return coefs, dependencies
        else:
            return coefs


def compute_micro_coefficients_old( conf, options, ret_all = False ):
    """Dependencies must be listed in a correct order."""
    aux = get_default_attr( conf, 'evp_equivalence', None )
    evp_equivalence = build_evp_equivalence( aux )

    opts = conf.options
    if hasattr( opts, 'post_process_hook' ) and opts.post_process_hook is not None:
        # User postprocessing.
        post_process_hook = getattr( conf.funmod, opts.post_process_hook )
    else:
        post_process_hook = None
    file_per_var = get_default_attr( opts, 'file_per_var', True )

    problem = ProblemDefinition.from_conf( conf,
                                           init_variables = False,
                                           init_equations = False )

    dependencies = {}
    cache_evp = {}
    
    coefs = pfdpm.Coefficients()
    coefs.filename = conf._filename
    
    coef_info = getattr( conf, opts.coef_info )
    for coef_name, cargs in coef_info.iteritems():
        output( 'computing %s...' % coef_name )
        requires = cargs.get( 'requires', [] )
        for req in requires:
            if dependencies.has_key( req ) and (dependencies[req] is not None):
                continue

            output( 'computing dependency %s...' % req )

            if req == 'pis':
                dependencies['pis'] = pfdpm.create_pis( problem,
                                                       conf.variables, 'uc' )

            elif req == 'corrs_rs':
                ##
                # Solve steady_rs.
                solve = pfdpm.solve_steady_correctors_rs
                aux = solve( problem, conf.ebcs, conf.epbcs, conf.variables,
                             conf.equations_steady_rs,
                             conf.ebc_sets['corrs_rs'], conf.epbc_sets['corrs_rs'],
                             dependencies['pis'], opts.file_conf[req],
                             post_process_hook, file_per_var )
                dependencies['corrs_rs'] = aux

            elif req == 'corrs_time_rs':
                evp = get_evp( req, cache_evp, problem, conf,
                              equivalence = evp_equivalence )

                dim = problem.domain.mesh.dim
                ts = problem.get_time_solver().ts
                corrs_rs = dependencies['corrs_rs']

                output( 'time-variant rs correctors via eigensolutions...' )
                filenames_rs = nm.zeros( (dim, dim), dtype = nm.object )
                solve = pfdpm.make_t_correctors_via_evp
                for ir in range( dim ):
                    for ic in range( dim ):
                        filename =( opts.file_conf[req] % (ir,ic)) + '.h5'
                        filenames_rs[ir,ic] = filename
                        solve( problem, evp, -corrs_rs.states_rs[ir,ic], ts,
                               filename, post_process_hook = post_process_hook,
                               file_per_var = file_per_var )
                output( '...done' )

                if opts.check.get( 'time_correctors', False ):
                    output( 'verifying correctors %s...' % req )
                    verify = pfdpm.verify_t_correctors
                    ok = True
                    for ir in range( dim ):
                        for ic in range( dim ):
                            oo = verify( problem, evp.ebc, evp.epbc,
                                         conf.equations_time_debug,
                                         -corrs_rs.states_rs[ir,ic],
                                         filenames_rs[ir,ic] )
                            ok = ok and oo
                    output( '...done, ok: %s' % ok )


                dependencies['corrs_time_rs'] = filenames_rs

            elif req in ['corrs_alpha1', 'corrs_alpha2']:
                ##
                # Solve steady_alpha.
                alpha = int( req[-1] )
                solve = pfdpm.solve_steady_correctors_alpha
                aux = solve( problem, conf.ebcs, conf.epbcs, conf.variables,
                             conf.equations_steady_alpha,
                             conf.ebc_sets[req],
                             conf.epbc_sets[req],
                             opts.file_conf[req],
                             alpha, post_process_hook, file_per_var )
                dependencies[req] = aux

            elif req in ['corrs_time_alpha1', 'corrs_time_alpha2']:
                ts = problem.get_time_solver().ts
                alpha = int( req[-1] )
                corrs_alpha = dependencies['corrs_alpha%d' % alpha]

                output( 'time-variant alpha correctors via eigensolutions...' )
                solve = pfdpm.make_t_correctors_via_evp

                evp = get_evp( req, cache_evp, problem, conf,
                              equivalence = evp_equivalence )
                filename = (opts.file_conf[req] % alpha) + '.h5'
                solve( problem, evp, corrs_alpha.state, ts, filename,
                       alpha = alpha, post_process_hook = post_process_hook,
                       file_per_var = file_per_var )
                output( '...done' )

                if opts.check.get( 'time_correctors', False ):
                    output( 'verifying correctors %s...' % req )
                    verify = pfdpm.verify_t_correctors
                    ok = verify( problem, evp.ebc, evp.epbc,
                                 conf.equations_time_debug,
                                 corrs_alpha.state,
                                 filename )
                    output( '...done, ok: %s' % ok )
                dependencies[req] = filename

            elif req == 'corrs_pressure':
                ##
                # Solve steady_pressure.
                solve = pfdpm.solve_steady_correctors_pressure
                aux = solve( problem, conf.ebcs, conf.epbcs, conf.variables,
                             conf.equations_steady_pressure,
                             conf.ebc_sets['corrs_pressure'],
                             conf.epbc_sets['corrs_pressure'],
                             opts.file_conf[req], post_process_hook, file_per_var )
                dependencies['corrs_pressure'] = aux

            elif req == 'corrs_time_pressure':
                evp = get_evp( req, cache_evp, problem, conf,
                              equivalence = evp_equivalence )

                ts = problem.get_time_solver().ts
                corrs_pressure = dependencies['corrs_pressure']

                output( 'time-variant pressure correctors via eigensolutions...' )
                solve = pfdpm.make_t_correctors_via_evp
                filename = opts.file_conf[req] + '.h5'
                solve( problem, evp, corrs_pressure.state, ts,
                       filename, post_process_hook = post_process_hook,
                       file_per_var = file_per_var )

                output( '...done' )

                if opts.check.get( 'time_correctors', False ):
                    output( 'verifying correctors %s...' % req )
                    ok = pfdpm.verify_t_correctors( problem, evp.ebc, evp.epbc,
                                                  conf.equations_time_debug,
                                                  corrs_pressure.state,
                                                  filename )
                    output( '...done, ok: %s' % ok )

                dependencies['corrs_time_pressure'] = filename

            else:
                print 'unknown dependency: %s' % req
                raise ValueError

            output( '...done' )
                
        if coef_name == 'times':
            val = problem.get_time_solver().ts.times

        elif coef_name in ['C1', 'C2', 'C']:
            ##
            # Permeabilities related to channels.
            val = pfdpm.solve_permeability( problem,
                                           conf.ebcs, conf.epbcs,
                                           conf.variables,
                                           conf.equations_coefs, cargs )

        elif coef_name == 'VF':
            ##
            # Volume fractions based on 'uc' geometry.
            aux, volume_all = pfdpm.volumes( problem, conf.variables,
                                             conf.equations_coefs, cargs )
            vf = {}
            vv = 0.0
            for key, val in aux.iteritems():
                setattr( coefs, 'volume%s' % key, val )
                if key != volume_all:
                    vf[key] =  nm.array( val / aux[volume_all],
                                         dtype = nm.float64 )
                    vv += val
            assert_( abs( vv - coefs.volumeY ) < 1e-14 )
            val = vf

        elif coef_name == 'E':
            val = pfdpm.coef_e( problem,
                               dependencies['corrs_rs'],
                               dependencies['pis'],
                               conf.variables, conf.equations_coefs, cargs )

        elif coef_name in ['GStar1', 'GStar2']:
            alpha = req[-1]
            val = pfdpm.coef_g_star( problem, dependencies['corrs_alpha' + alpha],
                                   conf.variables, conf.equations_coefs, cargs )

        elif coef_name in ['GBar1', 'GBar2']:
            val = pfdpm.coef_g_bar( problem,
                                  conf.ebcs, conf.epbcs,
                                  conf.variables, conf.equations_coefs, cargs )

        elif coef_name in ['GPlus1', 'GPlus2']:
            alpha = req[-1]
            val = pfdpm.coef_g_plus( problem,
                                   dependencies['corrs_time_alpha' + alpha],
                                   conf.ebcs, conf.epbcs,
                                   conf.variables, conf.equations_coefs, cargs )

        elif coef_name == 'BiotE':
            val = pfdpm.coef_biot_e( problem,
                                   dependencies['corrs_rs'],
                                   conf.variables, conf.equations_coefs, cargs )

        elif coef_name == 'BiotMIR':
            val = pfdpm.coef_biot_mir( problem,
                                     dependencies['corrs_pressure'],
                                     conf.variables, conf.equations_coefs,
                                     cargs )
        elif coef_name == 'H':
            val = pfdpm.coef_h( problem,
                               dependencies['corrs_rs'],
                               dependencies['corrs_time_rs'],
                               conf.variables, conf.equations_coefs, cargs )

        elif coef_name == 'BiotH':
            val = pfdpm.coef_biot_h( problem,
                                   dependencies['corrs_time_rs'],
                                   conf.variables, conf.equations_coefs, cargs )

        elif coef_name == 'BiotH2':
            val = pfdpm.coef_biot_h2( problem,
                                    dependencies['corrs_time_rs'],
                                    dependencies['corrs_time_pressure'],
                                    dependencies['corrs_rs'],
                                    dependencies['corrs_pressure'],
                                    conf.variables, conf.equations_coefs, cargs )

        elif coef_name == 'BiotMHR':
            val = pfdpm.coef_biot_mhr( problem,
                                     dependencies['corrs_time_pressure'],
                                     conf.variables, conf.equations_coefs,
                                     cargs )

        elif coef_name in ['P1', 'P2']:
            alpha = req[-1]
            val = pfdpm.coef_p( problem,
                               dependencies['corrs_alpha' + alpha],
                               dependencies['pis'],
                               conf.variables, conf.equations_coefs, cargs )

        elif coef_name in ['R1', 'R2']:
            alpha = req[-1]
            val = pfdpm.coef_r( problem,
                               dependencies['corrs_time_alpha' + alpha],
                               dependencies['pis'],
                               conf.variables, conf.equations_coefs, cargs )

        else:
            print 'unknown coefficient: %s' % coef_name
            raise ValueError

        setattr( coefs, coef_name, val )
        output( '...done' )

    prec = nm.get_printoptions()[ 'precision']
    if hasattr( opts, 'print_digits' ):
        nm.set_printoptions( precision = opts.print_digits )
    print coefs
    nm.set_printoptions( precision = prec )
##     pause()

    if ret_all:
        return coefs, dependencies
    else:
        return coefs

##
# c: 13.06.2008, r: 22.06.2008
def verify_steady_solution( conf, options ):
    if not hasattr( conf, 'equations_steady' ):
        output( 'set "equations_steady" in the input!' )
        return False

    ok = True
    pb = ProblemDefinition.from_conf( conf )
    opts = conf.options
    if hasattr( opts, 'post_process_hook' ) and opts.post_process_hook is not None:
        # User postprocessing.
        pph = getattr( conf.funmod, opts.post_process_hook )
    else:
        pph = None
    state_t, data_t = solve_evolutionary_op( pb, options, post_process_hook = pph )

    pb.set_equations( conf.equations_steady )
    tsc = pb.ts_conf
    ts = pb.get_default_ts( tsc.t0, tsc.t1, tsc.dt, tsc.n_step, tsc.n_step - 1 )
    # To reassemble matrix with new equations.
    pb.set_linear( False )
    state_s, data_s = solve_stationary_op( pb, options, ts = ts,
                                       post_process_hook = pph )

    err = nla.norm( state_s - state_t )
    eps = get_default_attr( conf.options, 'steady_state_error', 1e-8 )
    if err > eps:
        ok = False

    output( 'error: %.2e' % err )

    aux = eval_term_op( state_s, 'dw_volume_wdot.i1.Omega( m.BiotMHI, q, p )',
                      pb, dw_mode = 'vector' )
    print nla.norm( aux )
    aux = eval_term_op( state_t, 'dw_volume_wdot.i1.Omega( m.BiotMHI, q, p )',
                      pb, dw_mode = 'vector' )
    print nla.norm( aux )

    import pylab
    pylab.plot( state_s )
    pylab.plot( state_t )
    pylab.show()

    return ok

usage = """%prog [options] filename_in"""

help = {
    'filename' :
    'basename of output file(s) [default: <basename of input file>]',
    'verify_steady' :
    'verify steady state',
}

##
# c: 11.07.2006, r: 18.06.2008
def main():
    version = open( op.join( init_sfepy.install_dir,
                             'VERSION' ) ).readlines()[0][:-1]

    parser = OptionParser( usage = usage, version = "%prog " + version )
    parser.add_option( "-o", "", metavar = 'filename',
                       action = "store", dest = "output_filename_trunk",
                       default = None, help = help['filename'] )
    parser.add_option( "-v", "--verify-steady",
                       action = "store_true", dest = "verify_steady",
                       default = False, help = help['verify_steady'] )
    (options, args) = parser.parse_args()

    if (len( args ) == 1):
        filename_in = args[0];
    else:
        parser.print_help(),
        return

    required, other = get_standard_keywords()


    if not options.verify_steady:
        required.remove( 'equations' )

    conf = ProblemConf.from_file( filename_in, required, other )
##     print conf
##     pause()

    if options.verify_steady:
        ok = verify_steady_solution( conf, options )
        if not ok:
            output( 'failed!' )
        else:
            output( 'ok!' )
    else:
        app = PorousMediaHomogenizationApp( conf, options, 'homogen:' )
        opts = conf.options
        if hasattr( opts, 'parametric_hook' ): # Parametric study.
            parametric_hook = getattr( conf, opts.parametric_hook )
            app.parametrize( parametric_hook )
        app()
#        coefs = compute_micro_coefficients_old( conf, options )

if __name__ == '__main__':
    main()
