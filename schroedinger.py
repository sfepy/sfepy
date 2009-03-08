#!/usr/bin/env python
"""
Electronic structure solver.

Type:

$ ./schroedinger.py

for usage and help.

"""
import os
import os.path as op
from optparse import OptionParser
from math import pi

from scipy.optimize import broyden3, bisection
from scipy.optimize.nonlin import excitingmixing

import init_sfepy
from sfepy.base.base import *
from sfepy.base.conf import ProblemConf, get_standard_keywords
from sfepy.base.la import eig, norm_l2_along_axis
from sfepy.base.log import Log
from sfepy.applications import SimpleApp
from sfepy.fem import eval_term_op, MeshIO, ProblemDefinition
import sfepy.base.ioutils as io
from sfepy.solvers import Solver

def guess_n_eigs( n_electron, n_eigs = None ):
    """
    Guess the number of eigenvalues (energies) to compute so that the smearing
    iteration converges. Passing n_eigs overrides the guess.
    """
    if n_eigs is not None: return n_eigs
    
    if n_electron > 2:
        n_eigs = int(1.2 * ((0.5 * n_electron) + 5))
    else:
        n_eigs = n_electron
    return n_eigs

def smear( energies, e_f, width, exponent ):
    energies = nm.atleast_1d( energies )

    e1, e2 = e_f - width, e_f + width

    val = nm.zeros_like( energies )

    ii = nm.where( energies <= e1 )[0]
    val[ii] = 2.0

    ii = nm.where( (energies > e1) & (energies <= e_f) )[0]
    val[ii] = 2.0 - nm.power((energies[ii] - e1) / width, exponent)

    ii = nm.where( (energies > e_f) & (energies < e2) )[0]
    val[ii] = 0.0 + nm.power((e2 - energies[ii]) / width, exponent)

    return val

def setup_smearing( eigs, n_electron, width = 0.1, exponent = 2.0 ):

    def objective( e_f ):
        r = nm.sum( smear( eigs, e_f, width, exponent ) ) - n_electron
#        print e_f, r
        return r

##     import pylab
##     x = nm.linspace(eigs[0], eigs[-1], 1000)
##     pylab.plot( x, smear( x, -3, width, exponent ) )
##     pylab.show()

##     import pylab
##     x = nm.linspace(eigs[0], eigs[-1], 1000)
##     pylab.plot( x, [objective(y) for y in x] )
##     pylab.show()

    try:
        e_f = bisection( objective, eigs[0], eigs[-1], xtol = 1e-12 )
    except AssertionError:
        e_f = None
##     print eigs
##     print e_f, e_f - width, e_f + width
##     print objective(e_f)
##     debug()

    def smear_tuned( energies ):
        return smear( energies, e_f, width, exponent )

##     import pylab
##     x = nm.linspace(eigs[0], eigs[-1], 1000)
##     pylab.plot( x, smear_tuned( x ) )
##     pylab.show()

    return e_f, smear_tuned

def update_state_to_output( out, pb, vec, name, fill_value = None ):
    """Convert 'vec' to output for saving and insert it into 'out'. """
    aux = pb.state_to_output( vec, fill_value )
    key = aux.keys()[0]
    out[name] = aux[key]

def wrap_function( function, args ):
    ncalls = [0]
    times = []
    results = []
    def function_wrapper( x ):
        ncalls[0] += 1
        tt = time.time()

        results[:] = function( x, *args )
        eigs, mtx_s_phi, vec_n, vec_vh, vec_vxc = results

        tt2 = time.time()
        if tt2 < tt:
            raise RuntimeError, '%f >= %f' % (tt, tt2)
        times.append( tt2 - tt )
        return vec_vh + vec_vxc - x
    return ncalls, times, function_wrapper, results

class SchroedingerApp( SimpleApp ):

    def process_options( options ):
        """Application options setup. Sets default values for missing
        non-compulsory options."""
        get = options.get_default_attr
        
        eigen_solver = get( 'eigen_solver', None,
                            'missing "eigensolver" in options!' )

        n_electron = get( 'n_electron', 5 )
        n_eigs = guess_n_eigs( n_electron, n_eigs = get( 'n_eigs', None ) )
        # None -> save all.
        save_eig_vectors = get( 'save_eig_vectors', None )

        log_filename = get( 'log_filename', 'log.txt' )
        iter_fig_name = get( 'iter_fig_name', 'iterations.pdf' )
        # Called after DFT iteration, can do anything, no return value.
        iter_hook = get( 'iter_hook', None )

        return Struct( **locals() )
    process_options = staticmethod( process_options )

    def process_dft_options( options ):
        """Application DFT options setup. Sets default values for missing
        non-compulsory options."""
        get = options.get_default_attr
        
        dft_solver = get( 'dft_solver', None,
                          'missing "dft" in options!' )

        return Struct( **locals() )
    process_dft_options = staticmethod( process_dft_options )

    def __init__( self, conf, options, output_prefix, **kwargs ):
        SimpleApp.__init__( self, conf, options, output_prefix,
                            init_equations = False )

        output_dir = self.problem.output_dir

        opts = self.app_options
        opts.log_filename = op.join( output_dir, opts.log_filename )
        opts.iter_fig_name = op.join( output_dir, opts.iter_fig_name )
        self.mesh_results_name = op.join( opts.output_dir,
                                          self.problem.get_output_name() )
        self.eig_results_name = op.join( opts.output_dir,
                                         self.problem.ofn_trunk + '_eigs.txt' )

    def setup_options( self ):
        SimpleApp.setup_options( self )
        opts = SchroedingerApp.process_options( self.conf.options )
        if self.options.dft:
            opts += SchroedingerApp.process_dft_options( self.conf.options )
        self.app_options += opts

        funmod = self.conf.funmod

        hook = self.app_options.iter_hook
        if hook is not None:
            hook = getattr( funmod, hook )
        self.iter_hook = hook
    
    def call( self ):
        options = self.options

        if options.dft:
            if options.plot:
                from sfepy.base.plotutils import pylab
                options.plot = pylab is not None

            evp = self.solve_eigen_problem_n()
        else:
            evp = self.solve_eigen_problem_1()

        output( "solution saved to %s" % self.problem.get_output_name() )
        output( "in %s" % self.app_options.output_dir )

	if self.post_process_hook_final is not None: # User postprocessing.
	    self.post_process_hook_final( self.problem, evp = evp )

        return evp

    def iterate( self, vec_vhxc, eig_solver, mtx_b, log, file_output,
                 n_electron = None ):
        from sfepy.physics import dft

        self.itercount += 1

        pb = self.problem
        opts = self.app_options

        n_electron = get_default( n_electron, opts.n_electron )
        
        pb.select_bcs( ebc_names = ['ZeroSurface'] )
        pb.update_materials( extra_mat_args = {'mat_v' : {'vhxc' : vec_vhxc}} )

        dummy = pb.create_state_vector()

        output( 'assembling lhs...' )
        tt = time.clock()
        mtx_a = eval_term_op( dummy, pb.conf.equations['lhs'], pb,
                              dw_mode = 'matrix', tangent_matrix = pb.mtx_a )
        output( '...done in %.2f s' % (time.clock() - tt) )

        assert_( nm.alltrue( nm.isfinite( mtx_a.data ) ) )

        output( 'computing the Ax=Blx Kohn-Sham problem...' )
        tt = time.clock()
        eigs, mtx_s_phi = eig_solver( mtx_a, mtx_b,
                                      opts.n_eigs, eigenvectors = True )
        output( '...done in %.2f s' % (time.clock() - tt) )
        n_eigs_ok = len(eigs)

        output( 'setting-up smearing...' )
        e_f, smear_tuned = setup_smearing( eigs, n_electron )
        output( 'Fermi energy:', e_f )
        if e_f is None:
            raise Exception("cannot find Fermi energy - exiting.")
        weights = smear_tuned(eigs)
        output( '...done' )
        
        if (weights[-1] > 1e-12):
            print n_eigs_ok
            print eigs
            raise Exception("not enough eigenvalues have converged - exiting.")

        output( "saving solutions, iter=%d..." % self.itercount )
        out = {}
        var_name = pb.variables.get_names( kind = 'state' )[0]
        for ii in xrange( n_eigs_ok ):
            vec_phi = pb.variables.make_full_vec( mtx_s_phi[:,ii] )
            update_state_to_output( out, pb, vec_phi, var_name+'%03d' % ii )
        name = op.join( opts.output_dir, "iter%d" % self.itercount )
        pb.save_state('.'.join((name, opts.output_format)), out=out)
        output( "...solutions saved" )

        vec_phi = nm.zeros_like( vec_vhxc )
        vec_n = nm.zeros_like( vec_vhxc )

        for ii in xrange( n_eigs_ok ):
            vec_phi = pb.variables.make_full_vec( mtx_s_phi[:,ii] )
            vec_n += weights[ii] * vec_phi ** 2

        vec_vxc = nm.zeros_like( vec_vhxc )
        for ii, val in enumerate( vec_n ):
            vec_vxc[ii] = dft.getvxc( val/(4*pi), 0 )

        pb.set_equations( pb.conf.equations_vh )
        pb.select_bcs( ebc_names = ['VHSurface'] )
        pb.variables['n'].data_from_data( vec_n )
        output( "solving Ax=b Poisson equation" )
        vec_vh = pb.solve()

        #sphere = eval_term_op( dummy, pb.conf.equations['sphere'], pb)
        #print sphere

        norm = nla.norm( vec_vh + vec_vxc )
        dnorm = abs(norm - self.norm_vhxc0)
        log( norm, max(dnorm,1e-20) ) # logplot of pure 0 fails.
        file_output( '%d: F(x) = |VH + VXC|: %f, abs(F(x) - F(x_prev)): %e'\
                     % (self.itercount, norm, dnorm) )

        file_output("-"*70)
        file_output('Fermi energy:', e_f)
        file_output("----------------------------------------")
        file_output(" #  |  eigs           | smearing")
        file_output("----|-----------------|-----------------")
        for ii in xrange( n_eigs_ok ):
            file_output("% 3d | %-15s | %-15s" % (ii+1, eigs[ii], weights[ii]))
        file_output("----------------------------------------")
        file_output("|N|:   ", nla.norm(vec_n))
        file_output("|V_H|: ", nla.norm(vec_vh))
        file_output("|V_XC|:", nla.norm(vec_vxc))
        file_output("-"*70)

	if self.iter_hook is not None: # User postprocessing.
            data = Struct( eigs = eigs, mtx_s_phi = mtx_s_phi,
                           vec_n = vec_n, vec_vh = vec_vh, vec_vxc = vec_vxc )
	    self.iter_hook( self.problem, data = data )

        self.norm_vhxc0 = norm
        
        return eigs, mtx_s_phi, vec_n, vec_vh, vec_vxc

    def solve_eigen_problem_n( self ):
        opts = self.app_options
        pb = self.problem

        dim = pb.domain.mesh.dim

        pb.set_equations( pb.conf.equations )
        pb.select_bcs( ebc_names = ['ZeroSurface'] )

        dummy = pb.create_state_vector()

        output( 'assembling rhs...' )
        tt = time.clock()
        mtx_b = eval_term_op( dummy, pb.conf.equations['rhs'], pb,
                              dw_mode = 'matrix',
                              tangent_matrix = pb.mtx_a.copy() )
        output( '...done in %.2f s' % (time.clock() - tt) )
        assert_( nm.alltrue( nm.isfinite( mtx_b.data ) ) )

        n_eigs = get_default( opts.n_eigs, pb.mtx_a.shape[0] )
##         mtx_a.save( 'a.txt', format='%d %d %.12f\n' )
##         mtx_b.save( 'b.txt', format='%d %d %.12f\n' )

        if self.options.plot:
            log_conf = {
                'is_plot' : True,
                'aggregate' : 1,
                'yscales' : ['linear', 'log'],
            }
        else:
            log_conf = {
                'is_plot' : False,
            }
        log =  Log.from_conf( log_conf, ([r'$|F(x)|$'], [r'$|F(x)-x|$']) )

        file_output = Output('').get_output_function( opts.log_filename,
                                                      combined = True )

        eig_conf = pb.get_solver_conf( opts.eigen_solver )
        eig_solver = Solver.any_from_conf( eig_conf )

        vec_vhxc = nm.zeros( (pb.variables.di.ptr[-1],), dtype = nm.float64 )

        self.norm_vhxc0 = nla.norm( vec_vhxc )
        self.itercount = 0
        aux = wrap_function( self.iterate,
                             (eig_solver, mtx_b, log, file_output) )
        ncalls, times, nonlin_v, results = aux

        # Create and call the DFT solver.
        dft_conf = pb.get_solver_conf( opts.dft_solver )
        dft_status = {}
        dft_solver = Solver.any_from_conf( dft_conf, fun = nonlin_v,
                                           status = dft_status )
        vec_vhxc = dft_solver( vec_vhxc )
        eigs, mtx_s_phi, vec_n, vec_vh, vec_vxc = results
        output( 'DFT iteration time [s]:', dft_status['time_stats'] )
        
        if self.options.plot:
            log( save_figure = opts.iter_fig_name, finished = True )
            pause()

        coor = pb.domain.get_mesh_coors()
        r2 = norm_l2_along_axis(coor, squared=True)
        vec_nr2 = vec_n * r2

        pb.select_bcs( ebc_names = ['ZeroSurface'] )
        mtx_phi = self.make_full( mtx_s_phi )
        out = {}
        update_state_to_output( out, pb, vec_n, 'n' )
        update_state_to_output( out, pb, vec_nr2, 'nr2' )
        update_state_to_output( out, pb, vec_vh, 'vh' )
        update_state_to_output( out, pb, vec_vxc, 'vxc' )
        self.save_results( eigs, mtx_phi, out = out )

        return Struct( pb = pb, eigs = eigs, mtx_phi = mtx_phi,
                       vec_n = vec_n, vec_nr2 = vec_nr2,
                       vec_vh = vec_vh, vec_vxc = vec_vxc )

    def solve_eigen_problem_1( self ):
        from sfepy.fem import Mesh

        options = self.options
        opts = self.app_options
        pb = self.problem

        dim = pb.domain.mesh.dim

        pb.set_equations( pb.conf.equations )
        pb.time_update()

        dummy = pb.create_state_vector()

        output( 'assembling lhs...' )
        tt = time.clock()
        mtx_a = eval_term_op( dummy, pb.conf.equations['lhs'], pb,
                              dw_mode = 'matrix',
                              tangent_matrix = pb.mtx_a )
        output( '...done in %.2f s' % (time.clock() - tt) )

        output( 'assembling rhs...' )
        tt = time.clock()
        mtx_b = eval_term_op( dummy, pb.conf.equations['rhs'], pb,
                              dw_mode = 'matrix',
                              tangent_matrix = pb.mtx_a.copy() )
        output( '...done in %.2f s' % (time.clock() - tt) )

        n_eigs = get_default( opts.n_eigs, mtx_a.shape[0] )
##         mtx_a.save( 'a.txt', format='%d %d %.12f\n' )
##         mtx_b.save( 'b.txt', format='%d %d %.12f\n' )

        output( 'computing resonance frequencies...' )
        eig = Solver.any_from_conf( pb.get_solver_conf( opts.eigen_solver ) )
        eigs, mtx_s_phi = eig( mtx_a, mtx_b, n_eigs )
        output( '...done' )

        bounding_box = Mesh.from_file(pb.conf.filename_mesh).get_bounding_box()
        # this assumes a box (3D), or a square (2D):
        a = bounding_box[1][0] - bounding_box[0][0]
        E_exact = None
        if options.hydrogen or options.boron:
            if options.hydrogen:
                Z = 1
            elif options.boron:
                Z = 5
            if options.dim == 2:
                E_exact = [-float(Z)**2/2/(n-0.5)**2/4
                           for n in [1]+[2]*3+[3]*5 + [4]*8 + [5]*15]
            elif options.dim == 3:
                E_exact = [-float(Z)**2/2/n**2 for n in [1]+[2]*2**2+[3]*3**2 ]
        if options.well:
            if options.dim == 2:
                E_exact = [pi**2/(2*a**2)*x
                           for x in [2, 5, 5, 8, 10, 10, 13, 13,
                                     17, 17, 18, 20, 20 ] ]
            elif options.dim == 3:
                E_exact = [pi**2/(2*a**2)*x
                           for x in [3, 6, 6, 6, 9, 9, 9, 11, 11,
                                     11, 12, 14, 14, 14, 14, 14,
                                     14, 17, 17, 17] ]
        if options.oscillator:
            if options.dim == 2:
                E_exact = [1] + [2]*2 + [3]*3 + [4]*4 + [5]*5 + [6]*6
            elif options.dim == 3:
                E_exact = [float(1)/2+x for x in [1]+[2]*3+[3]*6+[4]*10 ]
        if E_exact is not None:
            print "a=%f" % a
            print "Energies:"
            print     "n      exact         FEM      error"

            for i, e in enumerate(eigs):
                from numpy import NaN
                if i < len(E_exact):
                    exact = E_exact[i]
                    err = 100*abs((exact - e)/exact)
                else:
                    exact = NaN
                    err = NaN
                print "%d:  %.8f   %.8f  %5.2f%%" % (i, exact, e, err)
        else:
            print eigs
##         import sfepy.base.plotutils as plu
##         plu.spy( mtx_b, eps = 1e-12 )
##         plu.pylab.show()
##         pause()

        mtx_phi = self.make_full( mtx_s_phi )
        self.save_results( eigs, mtx_phi )

        return Struct( pb = pb, eigs = eigs, mtx_phi = mtx_phi )

    def make_full( self, mtx_s_phi ):
        pb = self.problem

        mtx_phi = nm.empty( (pb.variables.di.ptr[-1], mtx_s_phi.shape[1]),
                            dtype = nm.float64 )
        for ii in xrange( mtx_s_phi.shape[1] ):
            mtx_phi[:,ii] = pb.variables.make_full_vec( mtx_s_phi[:,ii] )

        return mtx_phi

    def save_results( self, eigs, mtx_phi, out = None ):
        pb = self.problem

        save = self.app_options.save_eig_vectors
        out = get_default( out, {} )
        for ii in xrange( eigs.shape[0] ):
            if save is not None:
                if (ii > save[0]) and (ii < (n_eigs - save[1])): continue
            aux = pb.state_to_output( mtx_phi[:,ii] )
            key = aux.keys()[0]
            out[key+'%03d' % ii] = aux[key]

        pb.save_state(self.mesh_results_name, out=out)

        fd = open(self.eig_results_name, 'w')
        eigs.tofile( fd, ' ' )
        fd.close()


usage = """%prog [options] filename_in

Solver for electronic structure problems. 

You need to create a mesh (optionally specify a dimension):

    $ ./schroedinger.py --mesh --2d

and then pick a problem to solve, some examples below (the dimension is
determined by the mesh that you created above):

    $ ./schroedinger.py --hydrogen
    $ ./schroedinger.py --well
    $ ./schroedinger.py --dft

and visualize the result:

    $ paraview --data=mesh.vtk

"""

help = {
    'filename' : 'basename of output file(s) [default: %default.vtk]',
    'well' : "solve infinite potential well (particle in a box) problem",
    'oscillator' : "solve spherically symmetric linear harmonic oscillator (1 electron) problem",
    'hydrogen' : "solve the hydrogen atom",
    'boron' : "solve the boron atom with 1 electron",
    "mesh": "creates a mesh",
    "dim": "Create a 2D mesh, instead of the default 3D",
    "dft": "Do a DFT calculation",
    "plot": "plot convergence of DFT iterations (with --dft)",
}

def main():
    version = open( op.join( init_sfepy.install_dir,
                             'VERSION' ) ).readlines()[0][:-1]

    parser = OptionParser( usage = usage, version = "%prog " + version )
    parser.add_option( "--mesh",
                       action = "store_true", dest = "mesh",
                       default = False, help = help['mesh'] )
    parser.add_option( "--2d",
                       action = "store_true", dest = "dim2",
                       default = False, help = help['dim'] )
    parser.add_option( "-o", "", metavar = 'filename',
                       action = "store", dest = "output_filename_trunk",
                       default = "mesh", help = help['filename'] )
    parser.add_option( "--oscillator",
                       action = "store_true", dest = "oscillator",
                       default = False, help = help['oscillator'] )
    parser.add_option( "--well",
                       action = "store_true", dest = "well",
                       default = False, help = help['well'] )
    parser.add_option( "--hydrogen",
                       action = "store_true", dest = "hydrogen",
                       default = False, help = help['hydrogen'] )
    parser.add_option( "--boron",
                       action = "store_true", dest = "boron",
                       default = False, help = help['boron'] )
    parser.add_option( "--dft",
                       action = "store_true", dest = "dft",
                       default = False, help = help['dft'] )
    parser.add_option( "-p", "--plot",
                       action = "store_true", dest = "plot",
                       default = False, help = help['plot'] )

    options, args = parser.parse_args()

    if len( args ) == 1:
        filename_in = args[0];
    elif len( args ) == 0:
        if options.oscillator:
            dim = MeshIO.any_from_filename("tmp/mesh.vtk").read_dimension()
            if dim == 2:
                filename_in = "input/quantum/oscillator2d.py"
            else:
                assert_( dim == 3 )
                filename_in = "input/quantum/oscillator3d.py"
            options.dim = dim
            print "Dimension:", dim
        elif options.well:
            dim = MeshIO.any_from_filename("tmp/mesh.vtk").read_dimension()
            if dim == 2:
                filename_in = "input/quantum/well2d.py"
            else:
                assert_( dim == 3 )
                filename_in = "input/quantum/well3d.py"
            options.dim = dim
            print "Dimension:", dim
        elif options.hydrogen:
            dim = MeshIO.any_from_filename("tmp/mesh.vtk").read_dimension()
            if dim == 2:
                filename_in = "input/quantum/hydrogen2d.py"
            else:
                assert_( dim == 3 )
                filename_in = "input/quantum/hydrogen3d.py"
            options.dim = dim
            print "Dimension:", dim
        elif options.boron:
            dim = MeshIO.any_from_filename("tmp/mesh.vtk").read_dimension()
            if dim == 2:
                filename_in = "input/quantum/boron2d.py"
            else:
                assert_( dim == 3 )
                filename_in = "input/quantum/boron3d.py"
            options.dim = dim
            print "Dimension:", dim
        elif options.mesh:
            try:
                os.makedirs("tmp")
            except OSError, e:
                if e.errno != 17: # [Errno 17] File exists
                    raise
            if options.dim2:
                print "Dimension: 2"
                os.system("cp database/quantum/square.geo tmp/mesh.geo")
                os.system("gmsh -2 tmp/mesh.geo -format mesh")
                os.system("script/mesh_to_vtk.py tmp/mesh.mesh tmp/mesh.vtk")
            else:
                print "Dimension: 3"
                import geom
                from sfepy.fem.mesh import Mesh
                try:
                    from site_cfg import tetgen_path
                except ImportError:
                    tetgen_path = '/usr/bin/tetgen'
                os.system("gmsh -0 database/box.geo -o tmp/x.geo")
                g = geom.read_gmsh("tmp/x.geo")
                g.printinfo()
                geom.write_tetgen(g, "tmp/t.poly")
                geom.runtetgen("tmp/t.poly", a=0.03, Q=1.0, quadratic=False,
                        tetgenpath=tetgen_path)
                m = Mesh.from_file("tmp/t.1.node")
                m.write("tmp/mesh.vtk", io="auto")
            print "Mesh written to tmp/mesh.vtk"
            return
        elif options.dft:
            dim = MeshIO.any_from_filename("tmp/mesh.vtk").read_dimension()
            if dim == 2:
                filename_in = "input/quantum/dft2d.py"
            else:
                assert_( dim == 3 )
                filename_in = "input/quantum/dft3d.py"
            print "Dimension:", dim
            options.dim = dim
        else:
            parser.print_help()
            return
    else:
        parser.print_help()
        return

    required, other = get_standard_keywords()
    conf = ProblemConf.from_file( filename_in, required, other )

    app = SchroedingerApp(conf, options, 'schroedinger:')
    opts = conf.options
    if hasattr( opts, 'parametric_hook' ): # Parametric study.
        parametric_hook = getattr( conf, opts.parametric_hook )
        app.parametrize( parametric_hook )
    app()

if __name__ == '__main__':
    main()
