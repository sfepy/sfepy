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

from scipy.optimize import broyden3
from scipy.optimize.nonlin import excitingmixing

import init_sfepy
from sfepy.base.base import *
from sfepy.base.conf import ProblemConf, get_standard_keywords
from sfepy.base.la import eig
from sfepy.fem import eval_term_op, MeshIO, ProblemDefinition
import sfepy.base.ioutils as io
from sfepy.homogenization.phono import process_options
from sfepy.solvers import Solver

##
# c: 22.02.2008, r: 22.02.2008
def update_state_to_output( out, pb, vec, name, fill_value = None ):
    aux = pb.state_to_output( vec, fill_value )
    key = aux.keys()[0]
    out[name] = aux[key]

##
# c: 22.02.2008, r: 22.02.2008
def wrap_function( function, args ):
    ncalls = [0]
    times = []
    def function_wrapper( x ):
        ncalls[0] += 1
        tt = time.time()
        out = function( x, *args )
        eigs, mtx_s_phi, vec_n, vec_vh, vec_vxc = out
        print "-"*70
        print "eigs",eigs
        print "V_H",vec_vh
        print "V_XC",vec_vxc
        print "-"*70
        tt2 = time.time()
        if tt2 < tt:
            raise RuntimeError, '%f >= %f' % (tt, tt2)
        times.append( tt2 - tt )
        return vec_vh + vec_vxc
    return ncalls, times, function_wrapper

itercount = 0
##
# c: 22.02.2008, r: 03.03.2008
def iterate( vec_vhxc, pb, conf, eig_solver, n_eigs, mtx_b, n_electron = 5 ):
    global itercount
    itercount += 1
    from sfepy.physics import dft

    pb.update_materials( extra_mat_args = {'mat_v' : {'vhxc' : vec_vhxc}} )

    dummy = pb.create_state_vector()

    output( 'assembling lhs...' )
    tt = time.clock()
    mtx_a = eval_term_op( dummy, conf.equations['lhs'], pb,
                       dw_mode = 'matrix', tangent_matrix = pb.mtx_a )
    output( '...done in %.2f s' % (time.clock() - tt) )


    print 'computing the Ax=Blx Kohn-Sham problem...'
    eigs, mtx_s_phi = eig_solver( mtx_a, mtx_b, conf.options.n_eigs )

    if len(eigs) < n_electron:
        print len(eigs)
        print eigs
        raise Exception("Not enough eigenvalues have converged. Exiting.")

    print "saving solutions, iter=%d" % itercount
    out = {}
    var_name = pb.variables.get_names( kind = 'state' )[0]
    for ii in xrange( len(eigs) ):
        vec_phi = pb.variables.make_full_vec( mtx_s_phi[:,ii] )
        update_state_to_output( out, pb, vec_phi, var_name+'%03d' % ii )
    pb.save_state( "tmp/iter%d.vtk" % itercount, out = out )
    print "solutions saved"

    vec_phi = nm.zeros_like( vec_vhxc )
    vec_n = nm.zeros_like( vec_vhxc )
    for ii in xrange( n_electron ):
        vec_phi = pb.variables.make_full_vec( mtx_s_phi[:,ii] )
        vec_n += vec_phi ** 2


    vec_vxc = nm.zeros_like( vec_vhxc )
    for ii, val in enumerate( vec_n ):
        vec_vxc[ii] = dft.getvxc( val/(4*pi), 0 )

    pb.set_equations( conf.equations_vh )
    pb.time_update()
    pb.variables['n'].data_from_data( vec_n )
    print "Solving Ax=b Poisson equation"
    vec_vh = pb.solve()

    #sphere = eval_term_op( dummy, conf.equations['sphere'], pb)
    #print sphere

    return eigs, mtx_s_phi, vec_n, vec_vh, vec_vxc

##
# c: 01.02.2008, r: 03.03.2008
def solve_eigen_problem_n( conf, options ):


    pb = ProblemDefinition.from_conf( conf )
    dim = pb.domain.mesh.dim

    pb.time_update()

    dummy = pb.create_state_vector()

    output( 'assembling rhs...' )
    tt = time.clock()
    mtx_b = eval_term_op( dummy, conf.equations['rhs'], pb,
                       dw_mode = 'matrix', tangent_matrix = pb.mtx_a.copy() )
    output( '...done in %.2f s' % (time.clock() - tt) )

    #mtxA.save( 'tmp/a.txt', format='%d %d %.12f\n' )
    #mtxB.save( 'tmp/b.txt', format='%d %d %.12f\n' )

    try:
        n_eigs = conf.options.n_eigs
    except AttributeError:
        n_eigs = mtx_a.shape[0]

    if n_eigs is None:
        n_eigs = mtx_a.shape[0]

##     mtx_a.save( 'a.txt', format='%d %d %.12f\n' )
##     mtx_b.save( 'b.txt', format='%d %d %.12f\n' )

    eig_conf = pb.get_solver_conf( conf.options.eigen_solver )
    eig_solver = Solver.any_from_conf( eig_conf )
    vec_vhxc = nm.zeros( (pb.variables.di.ptr[-1],), dtype = nm.float64 )
    ncalls, times, nonlin_v = wrap_function( iterate,
                                           (pb, conf, eig_solver, n_eigs, mtx_b) )

    vec_vhxc = broyden3( nonlin_v, vec_vhxc, verbose = True )
    out = iterate( vec_vhxc, pb, conf, eig_solver, n_eigs, mtx_b )
    eigs, mtx_s_phi, vec_n, vec_vh, vec_vxc = out

    coor = pb.domain.get_mesh_coors()
    r = coor[:,0]**2 + coor[:,1]**2 + coor[:,2]**2
    vec_nr2 = vec_n * r

    n_eigs = eigs.shape[0]
    opts = process_options( conf.options, n_eigs )

    mtx_phi = nm.empty( (pb.variables.di.ptr[-1], mtx_s_phi.shape[1]),
                       dtype = nm.float64 )
    for ii in xrange( n_eigs ):
        mtx_phi[:,ii] = pb.variables.make_full_vec( mtx_s_phi[:,ii] )

    out = {}
    for ii in xrange( n_eigs ):
        if opts.save is not None:
            if (ii > opts.save[0]) and (ii < (n_eigs - opts.save[1])): continue
        aux = pb.state_to_output( mtx_phi[:,ii] )
        key = aux.keys()[0]
        out[key+'%03d' % ii] = aux[key]

    update_state_to_output( out, pb, vec_n, 'n' )
    update_state_to_output( out, pb, vec_nr2, 'nr2' )
    update_state_to_output( out, pb, vec_vh, 'vh' )
    update_state_to_output( out, pb, vec_vxc, 'vxc' )

    ofn_trunk = options.output_filename_trunk
    pb.domain.mesh.write( ofn_trunk + '.vtk', io = 'auto', out = out )

    fd = open( ofn_trunk + '_eigs.txt', 'w' )
    eigs.tofile( fd, ' ' )
    fd.close()

    return Struct( pb = pb, eigs = eigs, mtx_phi = mtx_phi )

##
# c: 01.02.2008, r: 03.03.2008
def solve_eigen_problem1( conf, options ):


    pb = ProblemDefinition.from_conf( conf )
    dim = pb.domain.mesh.dim

    pb.time_update()

    dummy = pb.create_state_vector()

    output( 'assembling lhs...' )
    tt = time.clock()
    mtx_a = eval_term_op( dummy, conf.equations['lhs'], pb,
                       dw_mode = 'matrix', tangent_matrix = pb.mtx_a )
    output( '...done in %.2f s' % (time.clock() - tt) )

    output( 'assembling rhs...' )
    tt = time.clock()
    mtx_b = eval_term_op( dummy, conf.equations['rhs'], pb,
                       dw_mode = 'matrix', tangent_matrix = pb.mtx_a.copy() )
    output( '...done in %.2f s' % (time.clock() - tt) )

    #mtxA.save( 'tmp/a.txt', format='%d %d %.12f\n' )
    #mtxB.save( 'tmp/b.txt', format='%d %d %.12f\n' )
    try:
        n_eigs = conf.options.n_eigs
    except AttributeError:
        n_eigs = mtx_a.shape[0]

    if n_eigs is None:
        n_eigs = mtx_a.shape[0]

##     mtx_a.save( 'a.txt', format='%d %d %.12f\n' )
##     mtx_b.save( 'b.txt', format='%d %d %.12f\n' )
    print 'computing resonance frequencies...'
    eig = Solver.any_from_conf( pb.get_solver_conf( conf.options.eigen_solver ) )
    eigs, mtx_s_phi = eig( mtx_a, mtx_b, conf.options.n_eigs )
    from sfepy.fem.mesh import Mesh
    bounding_box = Mesh.from_file("tmp/mesh.vtk").get_bounding_box()
    # this assumes a box (3D), or a square (2D):
    a = bounding_box[1][0] - bounding_box[0][0]
    E_exact = None
    if options.hydrogen or options.boron:
        if options.hydrogen:
            Z = 1
        elif options.boron:
            Z = 5
        if options.dim == 2:
            E_exact = [-float(Z)**2/2/(n-0.5)**2/4 for n in [1]+[2]*3+[3]*5 +\
                    [4]*8 + [5]*15]
        elif options.dim == 3:
            E_exact = [-float(Z)**2/2/n**2 for n in [1]+[2]*2**2+[3]*3**2 ]
    if options.well:
        if options.dim == 2:
            E_exact = [pi**2/(2*a**2)*x for x in [2, 5, 5, 8, 10, 10, 13, 13,
                17, 17, 18, 20, 20 ] ]
        elif options.dim == 3:
            E_exact = [pi**2/(2*a**2)*x for x in [3, 6, 6, 6, 9, 9, 9, 11, 11,
                11, 12, 14, 14, 14, 14, 14, 14, 17, 17, 17] ]
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
##     import sfepy.base.plotutils as plu
##     plu.spy( mtx_b, eps = 1e-12 )
##     plu.pylab.show()
##     pause()
    n_eigs = eigs.shape[0]
    opts = process_options( conf.options, n_eigs )

    mtx_phi = nm.empty( (pb.variables.di.ptr[-1], mtx_s_phi.shape[1]),
                       dtype = nm.float64 )
    for ii in xrange( n_eigs ):
        mtx_phi[:,ii] = pb.variables.make_full_vec( mtx_s_phi[:,ii] )

    out = {}
    for ii in xrange( n_eigs ):
        if opts.save is not None:
            if (ii > opts.save[0]) and (ii < (n_eigs - opts.save[1])): continue
        aux = pb.state_to_output( mtx_phi[:,ii] )
        key = aux.keys()[0]
        out[key+'%03d' % ii] = aux[key]

    ofn_trunk = options.output_filename_trunk
    pb.domain.mesh.write( ofn_trunk + '.vtk', io = 'auto', out = out )

    fd = open( ofn_trunk + '_eigs.txt', 'w' )
    eigs.tofile( fd, ' ' )
    fd.close()

    return Struct( pb = pb, eigs = eigs, mtx_phi = mtx_phi )


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
}

##
# c: 01.02.2008, r: 21.03.2008
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

    if options.dft:
        evp = solve_eigen_problem_n( conf, options )
    else:
        evp = solve_eigen_problem1( conf, options )

    print "Solution saved to %s.vtk" % options.output_filename_trunk

if __name__ == '__main__':
    main()
