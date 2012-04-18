import os
import time
from math import pi

import numpy as nm

from sfepy.base.base import Struct, output, get_default
from sfepy.applications import PDESolverApp
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

class SchroedingerApp(PDESolverApp):
    """
    Base application for electronic structure calculations.

    Subclasses should typically override `solve_eigen_problem()` method.

    This class allows solving only simple single electron problems,
    e.g. well, oscillator, hydrogen atom and boron atom with 1 electron.
    """

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
        PDESolverApp.__init__(self, conf, options, output_prefix,
                              init_equations=False)

    def setup_options(self):
        PDESolverApp.setup_options(self)
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
        self.mesh_results_name = os.path.join(opts.output_dir,
                                              self.problem.get_output_name())
        self.eig_results_name = os.path.join(opts.output_dir,
                                             self.problem.ofn_trunk
                                             + '_eigs.txt')

    def call(self):
        # This cannot be in __init__(), as parametric calls may change
        # the output directory.
        self.setup_output()

        evp = self.solve_eigen_problem()

        output("solution saved to %s" % self.problem.get_output_name())
        output("in %s" % self.app_options.output_dir)

        if self.post_process_hook_final is not None: # User postprocessing.
            self.post_process_hook_final(self.problem, evp=evp)

        return evp

    def solve_eigen_problem(self):
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

    def save_results(self, eigs, mtx_phi, out=None,
                     mesh_results_name=None, eig_results_name=None):
        mesh_results_name = get_default(mesh_results_name,
                                        self.mesh_results_name)
        eig_results_name = get_default(eig_results_name,
                                       self.eig_results_name)
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

        pb.save_state(mesh_results_name, out=out)

        fd = open(eig_results_name, 'w')
        eigs.tofile(fd, ' ')
        fd.close()
