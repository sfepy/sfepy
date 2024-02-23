"""
Eigenvalue problem solver application.
"""
from __future__ import absolute_import
import os

import numpy as nm

from sfepy.base.base import Struct, output, get_default
from sfepy.base.timing import Timer
from sfepy.applications import PDESolverApp
from sfepy.solvers import Solver
from sfepy.discrete.fem.meshio import convert_complex_output
from six.moves import range

class EVPSolverApp(PDESolverApp):
    """
    Solve an eigenvalue problem.
    """

    @staticmethod
    def process_options(options):
        """
        Application options setup. Sets default values for missing
        non-compulsory options.
        """
        get = options.get

        n_eigs = get('n_eigs', 5)
        eigs_only = get('eigs_only', False)
        return Struct(evps=get('evps', None, 'missing "evps" in options!'),
                      n_eigs=n_eigs,
                      eigs_only=eigs_only)

    def __init__(self, conf, options, output_prefix, **kwargs):
        PDESolverApp.__init__(self, conf, options, output_prefix,
                              init_equations=False)

    def setup_options(self):
        PDESolverApp.setup_options(self)
        opts = EVPSolverApp.process_options(self.conf.options)

        self.app_options += opts

    def setup_output(self):
        """
        Setup various file names for the output directory given by
        `self.problem.output_dir`.
        """
        output_dir = self.problem.output_dir

        opts = self.app_options
        opts.output_dir = output_dir
        self.mesh_results_name = self.problem.get_output_name()
        self.eig_results_name = os.path.join(opts.output_dir,
                                             self.problem.ofn_trunk
                                             + '_eigs.txt')

    def call(self, status=None):
        # This cannot be in __init__(), as parametric calls may change
        # the output directory.
        self.setup_output()

        evp = self.solve_eigen_problem()
        if status is not None:
            n_eigs = get_default(self.app_options.n_eigs, evp.vecs.shape[0])
            status.nls_status = Struct(conditions=n_eigs - len(evp.eigs))

        if self.post_process_hook_final is not None: # User postprocessing.
            self.post_process_hook_final(self.problem, evp=evp)

        return evp.pb, evp

    def solve_eigen_problem(self):
        opts = self.app_options
        pb = self.problem

        pb.set_equations(pb.conf.equations)
        pb.time_update()

        output('assembling lhs...')
        timer = Timer(start=True)
        mtx_a = pb.evaluate(pb.conf.equations['lhs'], mode='weak',
                            auto_init=True, dw_mode='matrix')
        output('...done in %.2f s' % timer.stop())

        if 'rhs' in pb.conf.equations:
            output('assembling rhs...')
            timer.start()
            mtx_b = pb.evaluate(pb.conf.equations['rhs'], mode='weak',
                                dw_mode='matrix')
            output('...done in %.2f s' % timer.stop())

        else:
            mtx_b = None

        _n_eigs = get_default(opts.n_eigs, mtx_a.shape[0])

        output('solving eigenvalue problem for {} values...'.format(_n_eigs))
        eig = Solver.any_from_conf(pb.get_solver_conf(opts.evps))
        if opts.eigs_only:
            eigs = eig(mtx_a, mtx_b, opts.n_eigs, eigenvectors=False)
            svecs = None

        else:
            eigs, svecs = eig(mtx_a, mtx_b, opts.n_eigs, eigenvectors=True)

        output('...done')

        vecs = self.make_full(svecs)
        self.save_results(eigs, vecs)

        return Struct(pb=pb, eigs=eigs, vecs=vecs)

    def make_full(self, svecs):
        if svecs is None: return None

        variables = self.problem.get_variables()

        vecs = nm.empty((variables.di.n_dof_total, svecs.shape[1]),
                        dtype=svecs.dtype)
        for ii in range(svecs.shape[1]):
            vecs[:,ii] = variables.make_full_vec(svecs[:,ii])

        return vecs

    def save_results(self, eigs, vecs, out=None,
                     mesh_results_name=None, eig_results_name=None):
        if vecs is not None:
            mesh_results_name = get_default(mesh_results_name,
                                            self.mesh_results_name)
            pb = self.problem
            pp_name = pb.conf.options.get('post_process_hook')
            pp = getattr(pb.conf.funmod,
                         pp_name if pp_name is not None else '',
                         lambda out, *args, **kwargs: out)

            out = get_default(out, {})
            variables = pb.set_default_state()
            for ii in range(eigs.shape[0]):
                variables.set_state(vecs[:, ii], force=True)
                aux = variables.create_output()
                aux2 = {}
                pp(aux2, pb, variables)
                aux.update(convert_complex_output(aux2))
                out.update({key + '%03d' % ii : aux[key] for key in aux})

            if aux.get('__mesh__') is not None:
                out['__mesh__'] = aux['__mesh__']

            pb.save_state(mesh_results_name, out=out)
            output('solution saved to %s' % mesh_results_name)

        eig_results_name = get_default(eig_results_name,
                                       self.eig_results_name)
        with open(eig_results_name, 'wb') as fd:
            if nm.iscomplexobj(eigs):
                nm.savetxt(fd, eigs, '% .18e % .18e')

            else:
                nm.savetxt(fd, eigs, '% .18e')

        output('eigenvalues saved to %s' % eig_results_name)
