#!/usr/bin/env python
"""
Electronic structure solver.

Type:

$ ./schroedinger.py

for usage and help.

"""
import os
import os.path as op
import time
from optparse import OptionParser
from math import pi

import numpy as nm

import sfepy
from sfepy.base.base import Struct, output, get_default
from sfepy.base.conf import ProblemConf, get_standard_keywords
from sfepy.base.ioutils import ensure_path
from sfepy.applications import SimpleApp
from sfepy.solvers import Solver

def guess_n_eigs(n_electron, n_eigs=None):
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

class SchroedingerApp(SimpleApp):

    @staticmethod
    def process_options(options):
        """
        Application options setup. Sets default values for missing
        non-compulsory options.

        Options:

        save_eig_vectors : (from_largest, from_smallest) or None
            If None, save all.
        """
        get = options.get_default_attr

        n_electron = get('n_electron', 5)
        n_eigs = guess_n_eigs(n_electron, n_eigs=get('n_eigs', None))

        return Struct(eigen_solver=get('eigen_solver', None,
                                       'missing "eigen_solver" in options!'),
                      n_electron=n_electron,
                      n_eigs=n_eigs,
                      save_eig_vectors=get('save_eig_vectors', None))

    def __init__(self, conf, options, output_prefix, **kwargs):
        SimpleApp.__init__(self, conf, options, output_prefix,
                           init_equations=False)

    def setup_options(self):
        SimpleApp.setup_options(self)
        opts = SchroedingerApp.process_options(self.conf.options)

        self.app_options += opts

    def setup_output(self):
        """
        Setup various file names for the output directory given by
        `self.problem.output_dir`.
        """
        output_dir = self.problem.output_dir

        opts = self.app_options
        opts.output_dir = output_dir
        self.mesh_results_name = op.join(opts.output_dir,
                                         self.problem.get_output_name())
        self.eig_results_name = op.join(opts.output_dir,
                                        self.problem.ofn_trunk + '_eigs.txt')

    def call(self):
        # This cannot be in __init__(), as parametric calls may change
        # the output directory.
        self.setup_output()

        evp = self.solve_eigen_problem_1()

        output("solution saved to %s" % self.problem.get_output_name())
        output("in %s" % self.app_options.output_dir)

        if self.post_process_hook_final is not None: # User postprocessing.
            self.post_process_hook_final(self.problem, evp=evp)

        return evp

    def solve_eigen_problem_1(self):
        options = self.options
        opts = self.app_options
        pb = self.problem

        dim = pb.domain.mesh.dim

        pb.set_equations(pb.conf.equations)
        pb.time_update()

        output('assembling lhs...')
        tt = time.clock()
        mtx_a = pb.evaluate(pb.conf.equations['lhs'], mode='weak',
                            auto_init=True, dw_mode='matrix')
        output('...done in %.2f s' % (time.clock() - tt))

        output('assembling rhs...')
        tt = time.clock()
        mtx_b = pb.evaluate(pb.conf.equations['rhs'], mode='weak',
                            dw_mode='matrix')
        output('...done in %.2f s' % (time.clock() - tt))

        n_eigs = get_default(opts.n_eigs, mtx_a.shape[0])

        output('computing resonance frequencies...')
        eig = Solver.any_from_conf(pb.get_solver_conf(opts.eigen_solver))
        eigs, mtx_s_phi = eig(mtx_a, mtx_b, n_eigs)
        output('...done')

        bounding_box = pb.domain.mesh.get_bounding_box()
        # this assumes a box (3D), or a square (2D):
        a = bounding_box[1][0] - bounding_box[0][0]
        E_exact = None
        if options.hydrogen or options.boron:
            if options.hydrogen:
                Z = 1
            elif options.boron:
                Z = 5
            if dim == 2:
                E_exact = [-float(Z)**2/2/(n-0.5)**2/4
                           for n in [1]+[2]*3+[3]*5 + [4]*8 + [5]*15]
            elif dim == 3:
                E_exact = [-float(Z)**2/2/n**2 for n in [1]+[2]*2**2+[3]*3**2 ]
        if options.well:
            if dim == 2:
                E_exact = [pi**2/(2*a**2)*x
                           for x in [2, 5, 5, 8, 10, 10, 13, 13,
                                     17, 17, 18, 20, 20 ] ]
            elif dim == 3:
                E_exact = [pi**2/(2*a**2)*x
                           for x in [3, 6, 6, 6, 9, 9, 9, 11, 11,
                                     11, 12, 14, 14, 14, 14, 14,
                                     14, 17, 17, 17] ]
        if options.oscillator:
            if dim == 2:
                E_exact = [1] + [2]*2 + [3]*3 + [4]*4 + [5]*5 + [6]*6
            elif dim == 3:
                E_exact = [float(1)/2+x for x in [1]+[2]*3+[3]*6+[4]*10 ]
        if E_exact is not None:
            output("a=%f" % a)
            output("Energies:")
            output("n      exact         FEM      error")

            for i, e in enumerate(eigs):
                from numpy import NaN
                if i < len(E_exact):
                    exact = E_exact[i]
                    err = 100*abs((exact - e)/exact)
                else:
                    exact = NaN
                    err = NaN
                output("%d:  %.8f   %.8f  %5.2f%%" % (i, exact, e, err))
        else:
            output(eigs)

        mtx_phi = self.make_full(mtx_s_phi)
        self.save_results(eigs, mtx_phi)

        return Struct(pb=pb, eigs=eigs, mtx_phi=mtx_phi)

    def make_full(self, mtx_s_phi):
        variables = self.problem.get_variables()

        mtx_phi = nm.empty((variables.di.ptr[-1], mtx_s_phi.shape[1]),
                           dtype=nm.float64)
        for ii in xrange(mtx_s_phi.shape[1]):
            mtx_phi[:,ii] = variables.make_full_vec(mtx_s_phi[:,ii])

        return mtx_phi

    def save_results(self, eigs, mtx_phi, out=None):
        pb = self.problem

        save = self.app_options.save_eig_vectors
        n_eigs = self.app_options.n_eigs
        out = get_default(out, {})
        state = pb.create_state()
        for ii in xrange(eigs.shape[0]):
            if save is not None:
                if (ii > save[0]) and (ii < (n_eigs - save[1])): continue
            state.set_full(mtx_phi[:,ii])
            aux = state.create_output_dict()
            key = aux.keys()[0]
            out[key+'%03d' % ii] = aux[key]

        pb.save_state(self.mesh_results_name, out=out)

        fd = open(self.eig_results_name, 'w')
        eigs.tofile(fd, ' ')
        fd.close()

def fix_path(filename):
    return op.join(sfepy.data_dir, filename)

usage = """%prog [options] [filename_in]

Solver for electronic structure problems. 

You need to create a mesh (optionally specify a dimension):

    $ ./schroedinger.py --create-mesh --2d

and then pick a problem to solve, some examples below (the dimension is
determined by the mesh that you created above):

    $ ./schroedinger.py --hydrogen
    $ ./schroedinger.py --well
    $ ./schroedinger.py --boron
    $ ./schroedinger.py --oscillator

and visualize the result:

- using Mayavi

  - 2D:

    $ ./postproc.py mesh.vtk

  - 3D:

    $ ./postproc.py mesh.vtk --3d

- using ParaView

    $ paraview --data=mesh.vtk
"""

help = {
    'conf' :
    'override problem description file items, written as python'
    ' dictionary without surrouding braces',
    'options' : 'override options item of problem description,'
    ' written as python dictionary without surrouding braces',
    'filename' :
    'basename of output file(s) [default: <basename of input file mesh>]',
    'well' :
    'solve infinite potential well (particle in a box) problem',
    'oscillator' :
    'solve spherically symmetric linear harmonic oscillator '
    '(1 electron) problem',
    'hydrogen' :
    'solve the hydrogen atom',
    'boron' :
    'solve the boron atom with 1 electron',
    'create_mesh' :
    'creates a mesh',
    'dim' :
    'Create a 2D mesh, instead of the default 3D',
    'mesh' :
    'override mesh file name. If given in form function(args), mesh is'
    ' generated on the fly.'
    ' Supported functions are cube(edge/2, nodes on edge) and'
    ' sphere(r, nodes density at start, end)',
    'mesh_dir' :
    'directory, where the mesh is generated by --mesh or --create-mesh'
    ' [default: %default]',
}

def main():
    parser = OptionParser(usage=usage, version='%prog ' + sfepy.__version__)
    parser.add_option('-c', '--conf', metavar='"key : value, ..."',
                      action='store', dest='conf', type='string',
                      default=None, help= help['conf'])
    parser.add_option('-O', '--options', metavar='"key : value, ..."',
                      action='store', dest='app_options', type='string',
                      default=None, help=help['options'])
    parser.add_option('-o', '', metavar='filename',
                      action='store', dest='output_filename_trunk',
                      default=None, help=help['filename'])
    parser.add_option('--create-mesh',
                      action='store_true', dest='create_mesh',
                      default=False, help=help['create_mesh'])
    parser.add_option('--2d',
                      action='store_true', dest='dim2',
                      default=False, help=help['dim'])
    parser.add_option('-m', '--mesh', metavar='filename',
                      action='store', dest='mesh',
                      default=None, help=help['mesh'])
    parser.add_option('--mesh-dir', metavar='dirname',
                      action='store', dest='mesh_dir',
                      default='tmp', help=help['mesh_dir'])
    parser.add_option('--oscillator',
                      action='store_true', dest='oscillator',
                      default=False, help=help['oscillator'])
    parser.add_option('--well',
                      action='store_true', dest='well',
                      default=False, help=help['well'])
    parser.add_option('--hydrogen',
                      action='store_true', dest='hydrogen',
                      default=False, help=help['hydrogen'])
    parser.add_option('--boron',
                      action='store_true', dest='boron',
                      default=False, help=help['boron'])

    options, args = parser.parse_args()

    if options.create_mesh and options.mesh:
        output('--create-mesh and --mesh options are mutually exclusive!')
        return

    if len(args) == 1:
        filename_in = args[0];
        auto_mesh_name = False

    elif len(args) == 0:
        auto_mesh_name = True

        mesh_filename = os.path.join(options.mesh_dir, 'mesh.vtk')
        ensure_path(mesh_filename)

        if options.oscillator:
            filename_in = fix_path("examples/quantum/oscillator.py")

        elif options.well:
            filename_in = fix_path("examples/quantum/well.py")

        elif options.hydrogen:
            filename_in = fix_path("examples/quantum/hydrogen.py")

        elif options.boron:
            filename_in = fix_path("examples/quantum/boron.py")

        elif options.create_mesh:
            output('generating mesh...')
            try:
                os.makedirs("tmp")
            except OSError, e:
                if e.errno != 17: # [Errno 17] File exists
                    raise
            if options.dim2:
                output("dimension: 2")
                gp = fix_path('meshes/quantum/square.geo')
                os.system("cp %s tmp/mesh.geo" % gp)
                os.system("gmsh -2 tmp/mesh.geo -format mesh")
                mtv = fix_path('script/mesh_to_vtk.py')
                os.system("%s tmp/mesh.mesh %s" % (mtv, mesh_filename))
            else:
                output("dimension: 3")
                import sfepy.geom as geom
                from sfepy.fem.mesh import Mesh
                try:
                    from site_cfg import tetgen_path
                except ImportError:
                    tetgen_path = '/usr/bin/tetgen'
                gp = fix_path('meshes/quantum/box.geo')
                os.system("gmsh -0 %s -o tmp/x.geo" % gp)
                g = geom.read_gmsh("tmp/x.geo")
                g.printinfo()
                geom.write_tetgen(g, "tmp/t.poly")
                geom.runtetgen("tmp/t.poly", a=0.03, Q=1.0,
                               quadratic=False, tetgenpath=tetgen_path)
                m = Mesh.from_file("tmp/t.1.node")
                m.write(mesh_filename, io="auto")
            output("...mesh written to %s" % mesh_filename)
            return

        else:
            parser.print_help()
            return

    else:
        parser.print_help()
        return

    required, other = get_standard_keywords()
    conf = ProblemConf.from_file_and_options(filename_in, options,
                                             required, other)

    if options.mesh:
        from sfepy.fem.mesh_generators import gen_mesh_from_string

        conf.filename_mesh = gen_mesh_from_string(options.mesh,
                                                  options.mesh_dir)

    elif auto_mesh_name and not sfepy.in_source_tree:
        conf.filename_mesh = mesh_filename
        conf.options.absolute_mesh_path = True

    app = SchroedingerApp(conf, options, 'schroedinger:')
    opts = conf.options
    if hasattr(opts, 'parametric_hook'): # Parametric study.
        parametric_hook = conf.get_function(opts.parametric_hook)
        app.parametrize(parametric_hook)
    app()

if __name__ == '__main__':
    main()
