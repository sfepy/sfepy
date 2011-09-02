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
import numpy.linalg as nla

try:
    from scipy.optimize import bisect
except ImportError:
    from scipy.optimize import bisection as bisect

import sfepy
from sfepy.base.base import Struct, Output, output, get_default, assert_, pause
from sfepy.base.conf import ProblemConf, get_standard_keywords
from sfepy.linalg import norm_l2_along_axis
from sfepy.base.log import Log
from sfepy.applications import SimpleApp
from sfepy.fem import Materials
from sfepy.fem.evaluate import eval_equations
from sfepy.fem.mappings import get_physical_qps, get_jacobian
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

def smear(energies, e_f, width, exponent):
    energies = nm.atleast_1d(energies)

    e1, e2 = e_f - width, e_f + width

    val = nm.zeros_like(energies)

    ii = nm.where(energies <= e1)[0]
    val[ii] = 2.0

    ii = nm.where((energies > e1) & (energies <= e_f))[0]
    val[ii] = 2.0 - nm.power((energies[ii] - e1) / width, exponent)

    ii = nm.where((energies > e_f) & (energies < e2))[0]
    val[ii] = 0.0 + nm.power((e2 - energies[ii]) / width, exponent)

    return val

def setup_smearing(eigs, n_electron, width=0.1, exponent=2.0):

    def objective(e_f):
        r = nm.sum(smear(eigs, e_f, width, exponent)) - n_electron

        return r

    try:
        e_f = bisect(objective, eigs[0], eigs[-1], xtol=1e-12)
    except AssertionError:
        e_f = None

    def smear_tuned(energies):
        return smear(energies, e_f, width, exponent)

    return e_f, smear_tuned

def update_state_to_output(out, pb, vec, name, fill_value=None):
    """Convert 'vec' to output for saving and insert it into 'out'. """
    state = pb.create_state()
    state.set_full(vec)
    aux = state.create_output_dict(fill_value=fill_value)
    key = aux.keys()[0]
    out[name] = aux[key]

def wrap_function(function, args, kwargs):
    ncalls = [0]
    times = []
    results = []
    def function_wrapper(x):
        ncalls[0] += 1
        tt = time.time()

        results[:] = function(x, *args, **kwargs)
        v_hxc_qp = results[-1]

        tt2 = time.time()
        if tt2 < tt:
            raise RuntimeError, '%f >= %f' % (tt, tt2)
        times.append(tt2 - tt)

        return v_hxc_qp.ravel() - x
    return ncalls, times, function_wrapper, results

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

    @staticmethod
    def process_dft_options(options):
        """
        Application DFT options setup. Sets default values for missing
        non-compulsory options.

        Options:

        save_dft_iterations : bool
            If True, save intermediate results during DFT iterations.
        init_hook : function
            The function called before DFT iterations, useful for
            defining additional keyword arguments passed to the
            eigenvalue problem solver.
        iter_hook : function
            The function called after each DFT iteration with no return
            value.
        iter_hook_final : function
            Like `iter_hook`, but called after the solver finishes.
        """
        get = options.get_default_attr

        return Struct(dft_solver=get('dft_solver', None,
                                     'missing "dft" in options!'),
                      log_filename=get('log_filename', 'log.txt'),
                      iter_fig_name=get('iter_fig_name', 'iterations.pdf'),
                      save_dft_iterations=get('save_dft_iterations', False),
                      init_hook=get('init_hook', None),
                      iter_hook=get('iter_hook', None),
                      iter_hook_final=get('iter_hook_final', None))

    def __init__(self, conf, options, output_prefix, **kwargs):
        SimpleApp.__init__(self, conf, options, output_prefix,
                           init_equations=False)

    def setup_options(self):
        SimpleApp.setup_options(self)
        opts = SchroedingerApp.process_options(self.conf.options)
        if self.options.dft:
            opts += SchroedingerApp.process_dft_options(self.conf.options)
            get = self.conf.get_function
            self.init_hook = get(opts.init_hook)
            self.iter_hook = get(opts.iter_hook)
            self.iter_hook_final = get(opts.iter_hook_final)

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
        if self.options.dft:
            opts.log_filename = op.join(output_dir, opts.log_filename)
            opts.iter_fig_name = op.join(output_dir, opts.iter_fig_name)

            # Restart file names for DFT.
            name = self.options.save_restart_filename
            if name == 'auto':
                name = op.join(output_dir,
                               self.problem.ofn_trunk + '_restart.h5')
            self.save_restart_filename = name

            name = self.options.load_restart_filename
            if name == 'auto':
                name = op.join(output_dir,
                               self.problem.ofn_trunk + '_restart.h5')
            self.load_restart_filename = name

    def call(self):
        options = self.options

        # This cannot be in __init__(), as parametric calls may change
        # the output directory.
        self.setup_output()

        if options.dft:
            if options.plot:
                from sfepy.base.plotutils import plt
                options.plot = plt is not None

            evp = self.solve_eigen_problem_n()
        else:
            evp = self.solve_eigen_problem_1()

        output("solution saved to %s" % self.problem.get_output_name())
        output("in %s" % self.app_options.output_dir)

        if self.post_process_hook_final is not None: # User postprocessing.
            self.post_process_hook_final(self.problem, evp=evp)

        return evp

    def _interp_to_nodes(self, v_qp):
        variable = self.problem.create_variables(['scalar'])['scalar']
        variable.data_from_qp(v_qp, self.problem.integrals['i1'])

        return variable()

    def iterate(self, v_hxc_qp, eig_solver,
                mtx_a_equations, mtx_a_variables, mtx_b, log, file_output,
                n_electron=None, **kwargs):
        from sfepy.physics import dft

        self.itercount += 1

        pb = self.problem
        opts = self.app_options

        n_electron = get_default(n_electron, opts.n_electron)

        sh = self.qp_shape

        v_hxc_qp = nm.array(v_hxc_qp, dtype=nm.float64)
        v_hxc_qp.shape = (sh[0] * sh[1],) + sh[2:]

        mat_v = Materials(mtx_a_equations.collect_materials())['mat_v']
        mat_v.set_extra_args(vhxc=v_hxc_qp)
        mat_v.time_update(None, mtx_a_equations, pb)

        v_hxc_qp.shape = sh

        v_ion_qp = mat_v.get_data(('Omega', 'i1'), 0, 'V_ion')
        for ig in xrange(1, pb.domain.shape.n_gr):
            _v_ion_qp = mat_v.get_data(('Omega', 'i1'), ig, 'V_ion')
            v_ion_qp = nm.concatenate((v_ion_qp, _v_ion_qp), axis=0)

        output('assembling lhs...')
        tt = time.clock()
        mtx_a = eval_equations(mtx_a_equations, mtx_a_variables,
                               mode='weak', dw_mode='matrix')
        output('...done in %.2f s' % (time.clock() - tt))

        assert_(nm.alltrue(nm.isfinite(mtx_a.data)))

        output('computing the Ax=Blx Kohn-Sham problem...')
        tt = time.clock()
        eigs, mtx_s_phi = eig_solver(mtx_a, mtx_b,
                                     opts.n_eigs, eigenvectors=True, **kwargs)
        output('...done in %.2f s' % (time.clock() - tt))
        n_eigs_ok = len(eigs)

        output('setting-up smearing...')
        e_f, smear_tuned = setup_smearing(eigs, n_electron)
        output('Fermi energy:', e_f)
        if e_f is None:
            raise Exception("cannot find Fermi energy - exiting.")
        weights = smear_tuned(eigs)
        output('...done')

        if (weights[-1] > 1e-12):
            output("last smearing weight is nonzero (%s eigs ok)!" % n_eigs_ok)

        if opts.save_dft_iterations:
            output("saving solutions, iter=%d..." % self.itercount)
            out = {}
            var_name = mtx_a_variables.get_names(kind='state')[0]
            for ii in xrange(n_eigs_ok):
                vec_phi = mtx_a_variables.make_full_vec(mtx_s_phi[:,ii])
                update_state_to_output(out, pb, vec_phi, var_name+'%03d' % ii)
            name = op.join(opts.output_dir, "iter%d" % self.itercount)
            pb.save_state('.'.join((name, opts.output_format)), out=out)
            output("...solutions saved")

        output('computing total charge...')
        tt = time.clock()
        aux = pb.create_evaluable('di_volume_integrate.i1.Omega(Psi)',
                                  mode='qp')
        psi_equations, psi_variables = aux
        var = psi_variables['Psi']

        n_qp = nm.zeros_like(v_hxc_qp)
        for ii in xrange(n_eigs_ok):
            vec_phi = mtx_a_variables.make_full_vec(mtx_s_phi[:,ii])
            var.data_from_any(vec_phi)

            phi_qp = eval_equations(psi_equations, psi_variables,
                                    mode='qp')
            n_qp += weights[ii] * (phi_qp ** 2)
        output('...done in %.2f s' % (time.clock() - tt))

        # Integrate charge density to get charge.
        det = get_jacobian(var.field, pb.integrals['i1'])
        charge = (det * n_qp).sum()

        vec_n = self._interp_to_nodes(n_qp)

        var.data_from_any(vec_n)
        charge_n = pb.evaluate('di_volume_integrate.i1.Omega(Psi)', Psi=var)

        ##
        # V_xc in quadrature points.
        v_xc_qp = nm.zeros((nm.prod(self.qp_shape),), dtype=nm.float64)
        for ii, val in enumerate(n_qp.flat):
            ## print ii, val
            v_xc_qp[ii] = dft.getvxc(val, 0, 0)
        assert_(nm.isfinite(v_xc_qp).all())
        v_xc_qp.shape = self.qp_shape

        mat_key = mat_v.datas.keys()[0]
        pb.set_equations(pb.conf.equations_vh)
        pb.select_bcs(ebc_names=['VHSurface'])
        pb.update_materials()

        qps = get_physical_qps(pb.domain.regions['Omega'], pb.integrals['i1'])
        mat_n_qp = {}
        for ig, ii in qps.indx.iteritems():
            mat_n_qp[ig] = {'N' : n_qp[ii]}

        output("solving Ax=b Poisson equation")
        mat_n = pb.get_materials()['mat_n']
        mat_n.reset()
        mat_n.set_all_data({mat_key : mat_n_qp})
        vec_v_h = pb.solve()()

        var.data_from_any(vec_v_h)
        v_h_qp = pb.evaluate('di_volume_integrate.i1.Omega(Psi)', Psi=var,
                             mode='qp')

        v_hxc_qp = v_h_qp + v_xc_qp
        norm = nla.norm(v_hxc_qp.ravel())
        dnorm = abs(norm - self.norm_v_hxc0)
        log(norm, max(dnorm,1e-20)) # logplot of pure 0 fails.
        file_output('%d: F(x) = |VH + VXC|: %f, abs(F(x) - F(x_prev)): %e'\
                    % (self.itercount, norm, dnorm))

        file_output("-"*70)
        file_output('Fermi energy:', e_f)
        file_output("----------------------------------------")
        file_output(" #  |  eigs           | smearing")
        file_output("----|-----------------|-----------------")
        for ii in xrange(n_eigs_ok):
            file_output("% 3d | %-15s | %-15s" % (ii+1, eigs[ii], weights[ii]))
        file_output("----------------------------------------")
        file_output("charge_qp: ", charge)
        file_output("charge_n:  ", charge_n)
        file_output("----------------------------------------")
        file_output("|N|:       ", nla.norm(n_qp.ravel()))
        file_output("|V_H|:     ", nla.norm(v_h_qp.ravel()))
        file_output("|V_XC|:    ", nla.norm(v_xc_qp.ravel()))
        file_output("|V_HXC|:   ", norm)

        if self.iter_hook is not None: # User postprocessing.
            pb.select_bcs(ebc_names=['ZeroSurface'])
            mtx_phi = self.make_full(mtx_s_phi)

            data = Struct(iteration=self.itercount,
                          eigs=eigs, weights=weights,
                          mtx_s_phi=mtx_s_phi, mtx_phi=mtx_phi,
                          vec_n=vec_n, vec_v_h=vec_v_h,
                          n_qp=n_qp, v_ion_qp=v_ion_qp, v_h_qp=v_h_qp,
                          v_xc_qp=v_xc_qp, file_output=file_output)
            self.iter_hook(self.problem, data=data)

        file_output("-"*70)

        self.norm_v_hxc0 = norm

        return eigs, weights, mtx_s_phi, vec_n, vec_v_h, \
               n_qp, v_ion_qp, v_h_qp, v_xc_qp, v_hxc_qp

    def solve_eigen_problem_n(self):
        opts = self.app_options
        pb = self.problem

        pb.set_equations(pb.conf.equations)
        pb.select_bcs(ebc_names=['ZeroSurface'])

        output('assembling rhs...')
        tt = time.clock()
        mtx_b = pb.evaluate(pb.conf.equations['rhs'], mode='weak',
                            auto_init=True, dw_mode='matrix')
        output('...done in %.2f s' % (time.clock() - tt))
        assert_(nm.alltrue(nm.isfinite(mtx_b.data)))

        ## mtx_b.save('b.txt', format='%d %d %.12f\n')

        aux = pb.create_evaluable(pb.conf.equations['lhs'], mode='weak',
                                  dw_mode='matrix')
        mtx_a_equations, mtx_a_variables = aux

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
        log =  Log.from_conf(log_conf, ([r'$|F(x)|$'], [r'$|F(x)-x|$']))

        file_output = Output('', opts.log_filename, combined=True)

        eig_conf = pb.get_solver_conf(opts.eigen_solver)
        eig_solver = Solver.any_from_conf(eig_conf)

        # Just to get the shape. Assumes one element group only!!!
        v_hxc_qp = pb.evaluate('di_volume_integrate.i1.Omega(Psi)',
                               mode='qp')
        self.qp_shape = v_hxc_qp.shape

        if self.load_restart_filename is not None:
            # Load a restart file.
            aux = self.load_dict(self.load_restart_filename)['V_hxc_qp']
            assert_(aux.shape == self.qp_shape)
            v_hxc_qp[:] = aux

        else:
            v_hxc_qp.fill(0.0)

        self.norm_v_hxc0 = nla.norm(v_hxc_qp)

        kwargs = self.init_hook(pb) if self.init_hook is not None else {}

        self.itercount = 0
        aux = wrap_function(self.iterate,
                            (eig_solver,
                             mtx_a_equations, mtx_a_variables,
                             mtx_b, log, file_output), kwargs)
        ncalls, times, nonlin_v, results = aux

        # Create and call the DFT solver.
        dft_conf = pb.get_solver_conf(opts.dft_solver)
        dft_status = {}
        dft_solver = Solver.any_from_conf(dft_conf,
                                          fun=nonlin_v,
                                          status=dft_status)
        v_hxc_qp = dft_solver(v_hxc_qp.ravel())

        v_hxc_qp = nm.array(v_hxc_qp, dtype=nm.float64)
        v_hxc_qp.shape = self.qp_shape
        eigs, weights, mtx_s_phi, vec_n, vec_v_h, \
              n_qp, v_ion_qp, v_h_qp, v_xc_qp, v_hxc_qp = results
        output('DFT iteration time [s]:', dft_status['time_stats'])

        fun = pb.functions['fun_v']
        variable = self.problem.create_variables(['scalar'])['scalar']
        field_coors = variable.field.get_coor()
        vec_v_ion = fun(None, field_coors,
                        mode='qp')['V_ion'].squeeze()

        vec_v_xc = self._interp_to_nodes(v_xc_qp)
        vec_v_hxc = self._interp_to_nodes(v_hxc_qp)
        vec_v_sum = self._interp_to_nodes(v_hxc_qp + v_ion_qp)

        r2 = norm_l2_along_axis(field_coors, squared=True)
        vec_nr2 = vec_n * r2

        pb.select_bcs(ebc_names=['ZeroSurface'])
        mtx_phi = self.make_full(mtx_s_phi)

        if self.save_restart_filename is not None:
            # Save a restart file.
            self.save_dict(self.save_restart_filename,
                           {'V_hxc_qp' : v_hxc_qp})

        if self.iter_hook_final is not None: # User postprocessing.
            data = Struct(iteration=self.itercount,
                          eigs=eigs, weights=weights,
                          mtx_s_phi=mtx_s_phi, mtx_phi=mtx_phi,
                          vec_n=vec_n, vec_v_h=vec_v_h,
                          n_qp=n_qp, v_ion_qp=v_ion_qp, v_h_qp=v_h_qp,
                          v_xc_qp=v_xc_qp, file_output=file_output)
            self.iter_hook_final(self.problem, data=data)

        out = {}
        update_state_to_output(out, pb, vec_n, 'n')
        update_state_to_output(out, pb, vec_nr2, 'nr2')
        update_state_to_output(out, pb, vec_v_h, 'V_h')
        update_state_to_output(out, pb, vec_v_xc, 'V_xc')
        update_state_to_output(out, pb, vec_v_ion, 'V_ion')
        update_state_to_output(out, pb, vec_v_hxc, 'V_hxc')
        update_state_to_output(out, pb, vec_v_sum, 'V_sum')
        self.save_results(eigs, mtx_phi, out=out)

        if self.options.plot:
            log(save_figure=opts.iter_fig_name)
            pause()
            log(finished=True)

        return Struct(pb=pb, eigs=eigs, mtx_phi=mtx_phi,
                      vec_n=vec_n, vec_nr2=vec_nr2,
                      vec_v_h=vec_v_h, vec_v_xc=vec_v_xc)

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

    $ ./schroedinger.py --mesh --2d

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
    'mesh' :
    'creates a mesh',
    'dim' :
    'Create a 2D mesh, instead of the default 3D',
    'dft' :
    'Do a DFT calculation (input file required)',
    'plot' :
    'plot convergence of DFT iterations (with --dft)',
    'save_restart' :
    'create restart file after DFT (with --dft)'
    ' special value "auto" is replaced by <output_dir>/<mesh_name>_restart.h5]',
    'load_restart' :
    'load restart file to initialize DFT (with --dft)'
    ' special value "auto" is replaced by <output_dir>/<mesh_name>_restart.h5]',
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
    parser.add_option('--mesh',
                      action='store_true', dest='mesh',
                      default=False, help=help['mesh'])
    parser.add_option('--2d',
                      action='store_true', dest='dim2',
                      default=False, help=help['dim'])
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
    parser.add_option('--dft',
                      action='store_true', dest='dft',
                      default=False, help=help['dft'])
    parser.add_option('-p', '--plot',
                      action='store_true', dest='plot',
                      default=False, help=help['plot'])
    parser.add_option('-r', '--save-restart', metavar='filename',
                      action='store', dest='save_restart_filename',
                      default=None, help=help['save_restart'])
    parser.add_option('-l', '--load-restart', metavar='filename',
                      action='store', dest='load_restart_filename',
                      default=None, help=help['load_restart'])

    options, args = parser.parse_args()

    if len(args) == 1:
        filename_in = args[0];
        auto_mesh_name = False

    elif len(args) == 0:
        auto_mesh_name = True

        if options.oscillator:
            filename_in = fix_path("examples/quantum/oscillator.py")

        elif options.well:
            filename_in = fix_path("examples/quantum/well.py")

        elif options.hydrogen:
            filename_in = fix_path("examples/quantum/hydrogen.py")

        elif options.boron:
            filename_in = fix_path("examples/quantum/boron.py")

        elif options.mesh:
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
                os.system("%s tmp/mesh.mesh tmp/mesh.vtk" % mtv)
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
                m.write("tmp/mesh.vtk", io="auto")
            output("...mesh written to tmp/mesh.vtk")
            return

        elif options.dft:
            output('the --dft option requires input file') 
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

    if auto_mesh_name and not sfepy.in_source_tree:
        conf.filename_mesh = "tmp/mesh.vtk"
        conf.options.absolute_mesh_path = True

    app = SchroedingerApp(conf, options, 'schroedinger:')
    opts = conf.options
    if hasattr(opts, 'parametric_hook'): # Parametric study.
        parametric_hook = conf.get_function(opts.parametric_hook)
        app.parametrize(parametric_hook)
    app()

if __name__ == '__main__':
    main()
