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
from sfepy.base.conf import ProblemConf, getStandardKeywords
from sfepy.base.la import eig
from sfepy.fem.evaluate import evalTermOP
from sfepy.fem.meshio import MeshIO
import sfepy.base.ioutils as io
from sfepy.fem.problemDef import ProblemDefinition
from sfepy.homogenization.phono import processOptions
from sfepy.solvers import Solver

##
# c: 22.02.2008, r: 22.02.2008
def updateStateToOutput( out, pb, vec, name, fillValue = None ):
    aux = pb.stateToOutput( vec, fillValue )
    key = aux.keys()[0]
    out[name] = aux[key]

##
# c: 22.02.2008, r: 22.02.2008
def wrapFunction( function, args ):
    ncalls = [0]
    times = []
    def function_wrapper( x ):
        ncalls[0] += 1
        tt = time.time()
        out = function( x, *args )
        eigs, mtxSPhi, vecN, vecVH, vecVXC = out
        print "-"*70
        print "eigs",eigs
        print "V_H",vecVH
        print "V_XC",vecVXC
        print "-"*70
        tt2 = time.time()
        if tt2 < tt:
            raise RuntimeError, '%f >= %f' % (tt, tt2)
        times.append( tt2 - tt )
        return vecVH + vecVXC
    return ncalls, times, function_wrapper

##
# c: 22.02.2008, r: 03.03.2008
def iterate( vecVHXC, pb, conf, eigSolver, nEigs, mtxB, nElectron = 5 ):
    from sfepy.physics import dft

    pb.updateMaterials( extraMatArgs = {'matV' : {'vhxc' : vecVHXC}} )

    dummy = pb.createStateVector()

    output( 'assembling lhs...' )
    tt = time.clock()
    mtxA = evalTermOP( dummy, conf.equations['lhs'], pb,
                       dwMode = 'matrix', tangentMatrix = pb.mtxA )
    output( '...done in %.2f s' % (time.clock() - tt) )


    print 'computing resonance frequencies...'
    eigs, mtxSPhi = eigSolver( mtxA, mtxB, conf.options.nEigs )

    if len(eigs) < nElectron:
        print len(eigs)
        print eigs
        raise Exception("Not enough eigenvalues have converged. Exitting.")

    vecPhi = nm.zeros_like( vecVHXC )
    vecN = nm.zeros_like( vecVHXC )
    for ii in xrange( nElectron ):
        vecPhi = pb.variables.makeFullVec( mtxSPhi[:,ii] )
        vecN += vecPhi ** 2


    vecVXC = nm.zeros_like( vecVHXC )
    for ii, val in enumerate( vecN ):
        vecVXC[ii] = dft.getvxc( val/(4*pi), 0 )

    pb.setEquations( conf.equations_vh )
    pb.timeUpdate()
    pb.variables['n'].dataFromData( vecN )
    vecVH = pb.solve()

    #sphere = evalTermOP( dummy, conf.equations['sphere'], pb)
    #print sphere

    return eigs, mtxSPhi, vecN, vecVH, vecVXC

##
# c: 01.02.2008, r: 03.03.2008
def solveEigenProblemN( conf, options ):


    pb = ProblemDefinition.fromConf( conf )
    dim = pb.domain.mesh.dim

    pb.timeUpdate()

    dummy = pb.createStateVector()

    output( 'assembling rhs...' )
    tt = time.clock()
    mtxB = evalTermOP( dummy, conf.equations['rhs'], pb,
                       dwMode = 'matrix', tangentMatrix = pb.mtxA.copy() )
    output( '...done in %.2f s' % (time.clock() - tt) )

    #mtxA.save( 'tmp/a.txt', format='%d %d %.12f\n' )
    #mtxB.save( 'tmp/b.txt', format='%d %d %.12f\n' )

    try:
        nEigs = conf.options.nEigs
    except AttributeError:
        nEigs = mtxA.shape[0]

    if nEigs is None:
        nEigs = mtxA.shape[0]

##     mtxA.save( 'a.txt', format='%d %d %.12f\n' )
##     mtxB.save( 'b.txt', format='%d %d %.12f\n' )

    eigConf = pb.getSolverConf( conf.options.eigenSolver )
    eigSolver = Solver.anyFromConf( eigConf )
    vecVHXC = nm.zeros( (pb.variables.di.ptr[-1],), dtype = nm.float64 )
    ncalls, times, nonlinV = wrapFunction( iterate,
                                           (pb, conf, eigSolver, nEigs, mtxB) )

    vecVHXC = broyden3( nonlinV, vecVHXC, verbose = True )
    out = iterate( vecVHXC, pb, conf, eigSolver, nEigs, mtxB )
    eigs, mtxSPhi, vecN, vecVH, vecVXC = out

    coor = pb.domain.getMeshCoors()
    r = coor[:,0]**2 + coor[:,1]**2 + coor[:,2]**2
    vecNr2 = vecN * r

    nEigs = eigs.shape[0]
    opts = processOptions( conf.options, nEigs )

    mtxPhi = nm.empty( (pb.variables.di.ptr[-1], mtxSPhi.shape[1]),
                       dtype = nm.float64 )
    for ii in xrange( nEigs ):
        mtxPhi[:,ii] = pb.variables.makeFullVec( mtxSPhi[:,ii] )

    out = {}
    for ii in xrange( nEigs ):
        if opts.save is not None:
            if (ii > opts.save[0]) and (ii < (nEigs - opts.save[1])): continue
        aux = pb.stateToOutput( mtxPhi[:,ii] )
        key = aux.keys()[0]
        out[key+'%03d' % ii] = aux[key]

    updateStateToOutput( out, pb, vecN, 'n' )
    updateStateToOutput( out, pb, vecNr2, 'nr2' )
    updateStateToOutput( out, pb, vecVH, 'vh' )
    updateStateToOutput( out, pb, vecVXC, 'vxc' )

    ofnTrunk = options.outputFileNameTrunk
    pb.domain.mesh.write( ofnTrunk + '.vtk', io = 'auto', out = out )

    fd = open( ofnTrunk + '_eigs.txt', 'w' )
    eigs.tofile( fd, ' ' )
    fd.close()

    return Struct( pb = pb, eigs = eigs, mtxPhi = mtxPhi )

##
# c: 01.02.2008, r: 03.03.2008
def solveEigenProblem1( conf, options ):


    pb = ProblemDefinition.fromConf( conf )
    dim = pb.domain.mesh.dim

    pb.timeUpdate()

    dummy = pb.createStateVector()

    output( 'assembling lhs...' )
    tt = time.clock()
    mtxA = evalTermOP( dummy, conf.equations['lhs'], pb,
                       dwMode = 'matrix', tangentMatrix = pb.mtxA )
    output( '...done in %.2f s' % (time.clock() - tt) )

    output( 'assembling rhs...' )
    tt = time.clock()
    mtxB = evalTermOP( dummy, conf.equations['rhs'], pb,
                       dwMode = 'matrix', tangentMatrix = pb.mtxA.copy() )
    output( '...done in %.2f s' % (time.clock() - tt) )

    #mtxA.save( 'tmp/a.txt', format='%d %d %.12f\n' )
    #mtxB.save( 'tmp/b.txt', format='%d %d %.12f\n' )
    try:
        nEigs = conf.options.nEigs
    except AttributeError:
        nEigs = mtxA.shape[0]

    if nEigs is None:
        nEigs = mtxA.shape[0]

##     mtxA.save( 'a.txt', format='%d %d %.12f\n' )
##     mtxB.save( 'b.txt', format='%d %d %.12f\n' )
    print 'computing resonance frequencies...'
    eig = Solver.anyFromConf( pb.getSolverConf( conf.options.eigenSolver ) )
    eigs, mtxSPhi = eig( mtxA, mtxB, conf.options.nEigs )
    print eigs
##     import sfepy.base.plotutils as plu
##     plu.spy( mtxB, eps = 1e-12 )
##     plu.pylab.show()
##     pause()
    nEigs = eigs.shape[0]
    opts = processOptions( conf.options, nEigs )

    mtxPhi = nm.empty( (pb.variables.di.ptr[-1], mtxSPhi.shape[1]),
                       dtype = nm.float64 )
    for ii in xrange( nEigs ):
        mtxPhi[:,ii] = pb.variables.makeFullVec( mtxSPhi[:,ii] )

    out = {}
    for ii in xrange( nEigs ):
        if opts.save is not None:
            if (ii > opts.save[0]) and (ii < (nEigs - opts.save[1])): continue
        aux = pb.stateToOutput( mtxPhi[:,ii] )
        key = aux.keys()[0]
        out[key+'%03d' % ii] = aux[key]

    ofnTrunk = options.outputFileNameTrunk
    pb.domain.mesh.write( ofnTrunk + '.vtk', io = 'auto', out = out )

    fd = open( ofnTrunk + '_eigs.txt', 'w' )
    eigs.tofile( fd, ' ' )
    fd.close()

    return Struct( pb = pb, eigs = eigs, mtxPhi = mtxPhi )


usage = """%prog [options] fileNameIn

Solver for electronic structure problems. 

You need to create a mesh (optionally specify a dimension):

    $ ./schroedinger --mesh -d2

and then pick a problem to solve, some examples below (specify the same
dimension as above):

    $ ./schroedinger --hydrogen -d2
    $ ./schroedinger --well
    $ ./schroedinger --dft

and visualize the result:

    $ paraview --data=mesh.vtk

"""

help = {
    'fileName' : 'basename of output file(s) [default: %default.vtk]',
    'well' : "solve infinite potential well (particle in a box) problem",
    'oscillator' : "solve spherically symmetric linear harmonic oscillator (1 electron) problem",
    'hydrogen' : "solve the hydrogen atom",
    "dim": "the dimensionality of the problem, either 2 or 3 [default: %default]",
    "mesh": "creates a mesh",
    "dft": "Uses a DFT solver"
}

##
# c: 01.02.2008, r: 21.03.2008
def main():
    version = open( op.join( init_sfepy.install_dir,
                             'VERSION' ) ).readlines()[0][:-1]

    parser = OptionParser( usage = usage, version = "%prog " + version )
    parser.add_option( "-o", "", metavar = 'fileName',
                       action = "store", dest = "outputFileNameTrunk",
                       default = "mesh", help = help['fileName'] )
    parser.add_option( "--oscillator",
                       action = "store_true", dest = "oscillator",
                       default = False, help = help['oscillator'] )
    parser.add_option( "--well",
                       action = "store_true", dest = "well",
                       default = False, help = help['well'] )
    parser.add_option( "-d", "--dim", type="int",
                       action = "store", dest = "dim",
                       default = 3, help = help['dim'] )
    parser.add_option( "--hydrogen",
                       action = "store_true", dest = "hydrogen",
                       default = False, help = help['hydrogen'] )
    parser.add_option( "--mesh",
                       action = "store_true", dest = "mesh",
                       default = False, help = help['mesh'] )
    parser.add_option( "--dft",
                       action = "store_true", dest = "dft",
                       default = False, help = help['dft'] )

    options, args = parser.parse_args()

    if len( args ) == 1:
        fileNameIn = args[0];
    elif len( args ) == 0:
        if options.oscillator:
            fileNameIn = "input/quantum/oscillator.py"
        elif options.well:
            fileNameIn = "input/quantum/well.py"
        elif options.hydrogen:
            # this is not reliable:
            #dim = MeshIO.anyFromFileName("tmp/mesh.vtk").read_dimension()
            if options.dim == 2:
                fileNameIn = "input/quantum/hydrogen2d.py"
            else:
                assert options.dim == 3
                fileNameIn = "input/quantum/hydrogen3d.py"
        elif options.mesh:
            try:
                os.makedirs("tmp")
            except OSError, e:
                if e.errno != 17: # [Errno 17] File exists
                    raise
            if options.dim == 2:
                os.system("cp database/square.geo tmp/mesh.geo")
                os.system("gmsh -2 tmp/mesh.geo -format mesh")
                os.system("script/mesh_to_vtk.py tmp/mesh.mesh tmp/mesh.vtk")
                print "Mesh written to tmp/mesh.vtk"
            else:
                assert options.dim == 3
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
                m = Mesh.fromFile("tmp/t.1.node")
                m.write("tmp/mesh.vtk", io="auto")
                print "Mesh written to tmp/mesh.vtk"
            return
        elif options.dft:
            if options.dim == 2:
                fileNameIn = "input/quantum/dft2d.py"
            else:
                assert options.dim == 3
                fileNameIn = "input/quantum/dft3d.py"
        else:
            parser.print_help()
            return
    else:
        parser.print_help()
        return

    required, other = getStandardKeywords()
    conf = ProblemConf.fromFile( fileNameIn, required, other )

    if options.dft:
        evp = solveEigenProblemN( conf, options )
    else:
        evp = solveEigenProblem1( conf, options )

    print "Solution saved to %s.vtk" % options.outputFileNameTrunk

if __name__ == '__main__':
    main()
