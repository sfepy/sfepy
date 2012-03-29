import os.path as op
import shutil
from copy import copy

import numpy as nm
import numpy.linalg as nla
import scipy as sc

from sfepy.base.base import output, set_defaults, get_default, assert_
from sfepy.base.base import Struct
from sfepy.base.conf import ProblemConf
from sfepy.fem import ProblemDefinition
from sfepy.fem.evaluate import assemble_by_blocks
from sfepy.homogenization.engine import HomogenizationEngine
from sfepy.homogenization.homogen_app import get_volume_from_options
from sfepy.homogenization.coefs_base import CoefDummy
from sfepy.applications import SimpleApp
from sfepy.solvers import Solver, eig
from sfepy.base.plotutils import plt

def make_save_hook(base_name, post_process_hook=None, file_per_var=None):
    def save_phono_correctors(state, problem, ir, ic):
        problem.save_state((base_name % (ir, ic)) + '.vtk', state,
                           post_process_hook=post_process_hook,
                           file_per_var=file_per_var)
    return save_phono_correctors

def try_set_defaults(obj, attr, defaults, recur=False):
    try:
        values = getattr(obj, attr)

    except:
        values = defaults

    else:
        if recur and isinstance(values, dict):
            for key, val in values.iteritems():
                set_defaults(val, defaults)

        else:
            set_defaults(values, defaults)

    return values

def transform_plot_data(datas, plot_transform, funmod):
    if plot_transform is not None:
        fun = getattr(funmod, plot_transform[0])

    dmin, dmax = 1e+10, -1e+10
    tdatas = []
    for data in datas:
        tdata = data.copy()
        if plot_transform is not None:
            tdata = fun(tdata, *plot_transform[1:])
        dmin = min(dmin, tdata.min())
        dmax = max(dmax, tdata.max())
        tdatas.append(tdata)
    dmin, dmax = min(dmax - 1e-8, dmin), max(dmin + 1e-8, dmax)
    return (dmin, dmax), tdatas

def plot_eigs(fig_num, plot_rsc, plot_labels, valid, freq_range, plot_range,
              show=False, clear=False, new_axes=False):
    """
    Plot resonance/eigen-frequencies.

    `valid` must correspond to `freq_range`

    resonances : red
    masked resonances: dotted red
    """
    if plt is None: return
    assert_(len(valid) == len(freq_range))

    fig = plt.figure(fig_num)
    if clear:
        fig.clf()
    if new_axes:
        ax = fig.add_subplot(111)
    else:
        ax = fig.gca()

    l0 = l1 = None
    for ii, f in enumerate(freq_range):
        if valid[ii]:
            l0 = ax.plot([f, f], plot_range, **plot_rsc['resonance'])[0]
        else:
            l1 = ax.plot([f, f], plot_range, **plot_rsc['masked'])[0]

    if l0:
        l0.set_label(plot_labels['resonance'])
    if l1:
        l1.set_label(plot_labels['masked'])

    if new_axes:
        ax.set_xlim([freq_range[0], freq_range[-1]])
        ax.set_ylim(plot_range)

    if show:
        plt.show()
    return fig

def plot_logs(fig_num, plot_rsc, plot_labels,
              freqs, logs, valid, freq_range, plot_range, squared,
              draw_eigs=True, show_legend=True, show=False,
              clear=False, new_axes=False):
    """
    Plot logs of min/middle/max eigs of a mass matrix.
    """
    if plt is None: return

    fig = plt.figure(fig_num)
    if clear:
        fig.clf()
    if new_axes:
        ax = fig.add_subplot(111)
    else:
        ax = fig.gca()

    if draw_eigs:
        plot_eigs(fig_num, plot_rsc, plot_labels, valid, freq_range,
                  plot_range)

    for ii, log in enumerate(logs):
        l1 = ax.plot(freqs[ii], log[:, -1], **plot_rsc['eig_max'])

        if log.shape[1] >= 2:
            l2 = ax.plot(freqs[ii], log[:, 0], **plot_rsc['eig_min'])
        else:
            l2 = None

        if log.shape[1] == 3:
            l3 = ax.plot(freqs[ii], log[:, 1], **plot_rsc['eig_mid'])
        else:
            l3 = None

    l1[0].set_label(plot_labels['eig_max'])
    if l2:
        l2[0].set_label(plot_labels['eig_min'])
    if l3:
        l3[0].set_label(plot_labels['eig_mid'])

    fmin, fmax = freqs[0][0], freqs[-1][-1]
    ax.plot([fmin, fmax], [0, 0], **plot_rsc['x_axis'])

    if squared:
        ax.set_xlabel(r'$\lambda$, $\omega^2$')
    else:
        ax.set_xlabel(r'$\sqrt{\lambda}$, $\omega$')

    ax.set_ylabel(plot_labels['y_axis'])

    if new_axes:
        ax.set_xlim([fmin, fmax])
        ax.set_ylim(plot_range)

    if show_legend:
        ax.legend()

    if show:
        plt.show()
    return fig

def plot_gap(ax, ii, f0, f1, kind, kind_desc, gmin, gmax, plot_range, plot_rsc):
    """
    Plot a single band gap as a rectangle.
    """
    def draw_rect(ax, x, y, rsc):
        ax.fill(nm.asarray(x)[[0,1,1,0]],
                 nm.asarray(y)[[0,0,1,1]],
                 **rsc)

    # Colors.
    strong = plot_rsc['strong_gap']
    weak = plot_rsc['weak_gap']
    propagation = plot_rsc['propagation']

    if kind == 'p':
        draw_rect(ax, (f0, f1), plot_range, propagation)
        info = [(f0, f1)]
    elif kind == 'w':
        draw_rect(ax, (f0, f1), plot_range, weak)
        info = [(f0, f1)]
    elif kind == 'wp':
        draw_rect(ax, (f0, gmin[1]), plot_range, weak)
        draw_rect(ax, (gmin[1], f1), plot_range, propagation)
        info = [(f0, gmin[1]), (gmin[1], f1)]
    elif kind == 's':
        draw_rect(ax, (f0, f1), plot_range, strong)
        info = [(f0, f1)]
    elif kind == 'sw':
        draw_rect(ax, (f0, gmax[1]), plot_range, strong)
        draw_rect(ax, (gmax[1], f1), plot_range, weak)
        info = [(f0, gmax[1]), (gmax[1], f1)]
    elif kind == 'swp':
        draw_rect(ax, (f0, gmax[1]), plot_range, strong)
        draw_rect(ax, (gmax[1], gmin[1]), plot_range, weak)
        draw_rect(ax, (gmin[1], f1), plot_range, propagation)
        info = [(f0, gmax[1]), (gmax[1], gmin[1]), (gmin[1], f1)]
    elif kind == 'is':
        draw_rect(ax, (gmin[1], gmax[1]), plot_range, strong)
        info = [(gmin[1], gmax[1])]
    elif kind == 'iw':
        draw_rect(ax, (gmin[1], gmax[1]), plot_range, weak)
        info = [(gmin[1], gmax[1])]
    else:
        output('impossible band gap combination:')
        output(gmin, gmax)
        raise ValueError

    output(ii, gmin[0], gmax[0], '%.8f' % f0, '%.8f' % f1)
    output(' -> %s\n    %s' %(kind_desc, info))

def plot_gaps(fig_num, plot_rsc, gaps, kinds, freq_range,
              plot_range, show=False, clear=False, new_axes=False):
    """
    Plot band gaps as rectangles.
    """
    if plt is None: return

    fig = plt.figure(fig_num)
    if clear:
        fig.clf()
    if new_axes:
        ax = fig.add_subplot(111)
    else:
        ax = fig.gca()

    for ii in xrange(len(freq_range) - 1):
        f0, f1 = freq_range[[ii, ii+1]]
        gap = gaps[ii]

        if isinstance(gap, list):
            for ig, (gmin, gmax) in enumerate(gap):
                kind, kind_desc = kinds[ii][ig]
                plot_gap(ax, ii, f0, f1, kind, kind_desc, gmin, gmax,
                         plot_range, plot_rsc)
        else:
            gmin, gmax = gap
            kind, kind_desc = kinds[ii]
            plot_gap(ax, ii, f0, f1, kind, kind_desc, gmin, gmax,
                     plot_range, plot_rsc)

    if new_axes:
        ax.set_xlim([freq_range[0], freq_range[-1]])
        ax.set_ylim(plot_range)

    if show:
        plt.show()
    return fig

def report_iw_cat(iw_dir, christoffel):
    output('incident wave direction:')
    output(iw_dir)
    output('Christoffel acoustic tensor:')
    output(christoffel)

class AcousticBandGapsApp(SimpleApp):
    """
    Application for computing acoustic band gaps.
    """

    @staticmethod
    def process_options(options):
        """
        Application options setup. Sets default values for missing
        non-compulsory options.
        """
        get = options.get_default_attr

        coefs_basic = get('coefs_basic', None,
                          'missing "coefs_basic" in options!')
        coefs_dispersion = get('coefs_dispersion', None,
                               'missing "coefs_dispersion" in options!')

        requirements = get('requirements', None,
                           'missing "requirements" in options!')

        default_plot_options = {'show' : True,'legend' : False,}

        aux = {
            'resonance' : 'eigenfrequencies',
            'masked' : 'masked eigenfrequencies',
            'eig_min' : 'min eig($M^*$)',
            'eig_mid' : 'mid eig($M^*$)',
            'eig_max' : 'max eig($M^*$)',
            'y_axis' : 'eigenvalues of mass matrix $M^*$',
        }
        plot_labels = try_set_defaults(options, 'plot_labels', aux, recur=True)

        aux = {
            'resonance' : 'eigenfrequencies',
            'masked' : 'masked eigenfrequencies',
            'eig_min' : r'$\kappa$(min)',
            'eig_mid' : r'$\kappa$(mid)',
            'eig_max' : r'$\kappa$(max)',
            'y_axis' : 'polarization angles',
        }
        plot_labels_angle = try_set_defaults(options, 'plot_labels_angle', aux)

        aux = {
            'resonance' : 'eigenfrequencies',
            'masked' : 'masked eigenfrequencies',
            'eig_min' : r'wave number (min)',
            'eig_mid' : r'wave number (mid)',
            'eig_max' : r'wave number (max)',
            'y_axis' : 'wave numbers',
        }
        plot_labels_wave = try_set_defaults(options, 'plot_labels_wave', aux)

        plot_rsc =  {
            'resonance' : {'linewidth' : 0.5, 'color' : 'r', 'linestyle' : '-' },
            'masked' : {'linewidth' : 0.5, 'color' : 'r', 'linestyle' : ':' },
            'x_axis' : {'linewidth' : 0.5, 'color' : 'k', 'linestyle' : '--' },
            'eig_min' : {'linewidth' : 0.5, 'color' : 'b', 'linestyle' : '--' },
            'eig_mid' : {'linewidth' : 0.5, 'color' : 'b', 'linestyle' : '-.' },
            'eig_max' : {'linewidth' : 0.5, 'color' : 'b', 'linestyle' : '-' },
            'strong_gap' : {'linewidth' : 0, 'facecolor' : (1, 1, 0.5) },
            'weak_gap' : {'linewidth' : 0, 'facecolor' : (1, 1, 1) },
            'propagation' : {'linewidth' : 0, 'facecolor' : (0.5, 1, 0.5) },
            'params' : {'axes.labelsize': 'large',
                        'text.fontsize': 'large',
                        'legend.fontsize': 'large',
                        'xtick.labelsize': 'large',
                        'ytick.labelsize': 'large',
                        'text.usetex': False},
        }
        plot_rsc = try_set_defaults(options, 'plot_rsc', plot_rsc)

        tensor_names = get('tensor_names', None,
                            'missing "tensor_names" in options!')

        volume = get('volume', None, 'missing "volume" in options!')

        return Struct(clear_cache=get('clear_cache', {}),
                      coefs_basic=coefs_basic,
                      coefs_dispersion=coefs_dispersion,
                      requirements=requirements,

                      incident_wave_dir=get('incident_wave_dir', None),
                      dispersion=get('dispersion', 'simple'),
                      dispersion_conf=get('dispersion_conf', None),
                      homogeneous=get('homogeneous', False),

                      eig_range=get('eig_range', None),
                      tensor_names=tensor_names,
                      volume=volume,

                      eig_vector_transform=get('eig_vector_transform', None),
                      plot_transform=get('plot_transform', None),
                      plot_transform_wave=get('plot_transform_wave', None),
                      plot_transform_angle=get('plot_transform_angle', None),

                      plot_options=get('plot_options', default_plot_options),

                      fig_name=get('fig_name', None),
                      fig_name_wave=get('fig_name_wave', None),
                      fig_name_angle=get('fig_name_angle', None),
                      fig_suffix=get('fig_suffix', '.pdf'),

                      plot_labels=plot_labels,
                      plot_labels_angle=plot_labels_angle,
                      plot_labels_wave=plot_labels_wave,
                      plot_rsc=plot_rsc)

    @staticmethod
    def process_options_pv(options):
        """
        Application options setup for phase velocity computation. Sets default
        values for missing non-compulsory options.
        """
        get = options.get_default_attr

        tensor_names = get('tensor_names', None,
                           'missing "tensor_names" in options!')

        volume = get('volume', None, 'missing "volume" in options!')

        return Struct(clear_cache=get('clear_cache', {}),

                      incident_wave_dir=get('incident_wave_dir', None),
                      dispersion=get('dispersion', 'simple'),
                      dispersion_conf=get('dispersion_conf', None),
                      homogeneous=get('homogeneous', False),
                      fig_suffix=get('fig_suffix', '.pdf'),

                      tensor_names=tensor_names,
                      volume=volume)

    def __init__(self, conf, options, output_prefix, **kwargs):
        SimpleApp.__init__(self, conf, options, output_prefix,
                           init_equations=False)

        self.setup_options()
        self.cached_coefs = None
        self.cached_iw_dir = None
        self.cached_christoffel = None
        self.cached_evp = None

        output_dir = self.problem.output_dir
        shutil.copyfile(conf._filename,
                        op.join(output_dir, op.basename(conf._filename)))

    def setup_options(self):
        SimpleApp.setup_options(self)

        if self.options.phase_velocity:
            process_options = AcousticBandGapsApp.process_options_pv
        else:
            process_options = AcousticBandGapsApp.process_options
        self.app_options += process_options(self.conf.options)

    def call(self):
        """
        In parametric runs, cached data (homogenized coefficients,
        Christoffel acoustic tensor and eigenvalue problem solution) are
        cleared according to 'clear_cache' aplication options.

        Example:

        clear_cache = {'cached_christoffel' : True, 'cached_evp' : True}
        """
        options = self.options

        for key, val in self.app_options.clear_cache.iteritems():
            if val and key.startswith('cached_'):
                setattr(self, key, None)

        if options.phase_velocity:
            # No band gaps in this case.
            return self.compute_phase_velocity()

        opts = self.app_options

        if options.detect_band_gaps:
            # Compute basic coefficients and data.
            coefs_name = opts.coefs_basic

        elif options.analyze_dispersion:
            # Compute basic + dispersion coefficients and data.
            conf = self.problem.conf

            coefs = copy(conf.get(opts.coefs_basic, {}))
            coefs.update(conf.get(opts.coefs_dispersion, {}))

            conf.coefs_all = coefs
            coefs_name = 'coefs_all'

        else:
            # Compute only the eigenvalue problems.
            conf = self.problem.conf

            names = [req for req in conf.get(opts.requirements, [''])
                     if req.startswith('evp')]
            coefs = {'dummy' : {'requires' : names,
                                'class' : CoefDummy,}}

            conf.coefs_dummy = coefs
            coefs_name = 'coefs_dummy'

        he_options = Struct(coefs=coefs_name, requirements=opts.requirements,
                            post_process_hook=self.post_process_hook)
        volume = get_volume_from_options(opts, self.problem)

        he = HomogenizationEngine(self.problem, options,
                                  app_options=he_options,
                                  volume=volume)
        coefs = he()

        if options.plot:
            if options.detect_band_gaps:
                self.plot_band_gaps(coefs)

            elif options.analyze_dispersion:
                self.plot_dispersion(coefs)

        return coefs

    def fix_eig_range(self, n_eigs):
        eig_range = get_default(self.app_options.eig_range, (0, n_eigs))
        if eig_range[-1] < 0:
            eig_range[-1] += n_eigs + 1

        assert_(eig_range[0] < (eig_range[1] - 1))
        assert_(eig_range[1] <= n_eigs)
        self.app_options.eig_range = eig_range

    def plot_band_gaps(self, coefs):
        opts = self.app_options

        bg_keys = [key for key in coefs.to_dict().keys()
                   if key.startswith('band_gaps')]

        plot_opts =  opts.plot_options
        plot_rsc = opts.plot_rsc

        plt.rcParams.update(plot_rsc['params'])

        for ii, key in enumerate(bg_keys):
            bg = coefs.get_default_attr(key)

            plot_labels =  opts.plot_labels.get(key, opts.plot_labels)

            plot_range, teigs = transform_plot_data(bg.logs.eigs,
                                                    opts.plot_transform,
                                                    self.conf.funmod)
            fig = plot_gaps(ii, plot_rsc, bg.gaps, bg.kinds,
                            bg.freq_range_margins, plot_range,
                            clear=True)
            fig = plot_logs(ii, plot_rsc, plot_labels, bg.logs.freqs, teigs,
                            bg.valid[bg.eig_range],
                            bg.freq_range_initial,
                            plot_range, False,
                            show_legend=plot_opts['legend'],
                            new_axes=True)

            if opts.fig_name is not None:
                bg_name = key.replace('band_gaps', '')
                if not bg_name.startswith('_'):
                    bg_name = '_' + bg_name

                fig_name = opts.fig_name + bg_name + opts.fig_suffix

                fig.savefig(op.join(self.problem.output_dir, fig_name))

        if plot_opts['show']:
            plt.show()

    def solve_eigen_problem(self, ofn_trunk=None, post_process_hook=None):

        if self.cached_evp is not None:
            return self.cached_evp

        problem = self.problem
        ofn_trunk = get_default(ofn_trunk, problem.ofn_trunk,
                                'output file name trunk missing!')
        post_process_hook = get_default(post_process_hook,
                                        self.post_process_hook)

        conf = self.conf

        eig_problem = self.app_options.eig_problem
        if eig_problem in ['simple', 'simple_liquid']:
            problem.set_equations(conf.equations)
            problem.time_update()
            problem.update_materials()

            mtx_a = problem.evaluate(conf.equations['lhs'], mode='weak',
                                     auto_init=True, dw_mode='matrix')

            mtx_m = problem.evaluate(conf.equations['rhs'], mode='weak',
                                     dw_mode='matrix')

        elif eig_problem == 'schur':
            # A = K + B^T D^{-1} B.
            mtx = assemble_by_blocks(conf.equations, self.problem,
                                     ebcs=conf.ebcs, epbcs=conf.epbcs)
            problem.set_equations(conf.equations)
            problem.time_update()
            problem.update_materials()

            ls = Solver.any_from_conf(problem.ls_conf,
                                      presolve=True, mtx=mtx['D'])

            mtx_b, mtx_m = mtx['B'], mtx['M']
            mtx_dib = nm.empty(mtx_b.shape, dtype=mtx_b.dtype)
            for ic in xrange(mtx_b.shape[1]):
                mtx_dib[:,ic] = ls(mtx_b[:,ic].toarray().squeeze())
            mtx_a = mtx['K'] + mtx_b.T * mtx_dib

        else:
            raise NotImplementedError

##     from sfepy.base.plotutils import spy, plt
##     spy(mtx_b, eps = 1e-12)
##     plt.show()
##     mtx_a.save('a.txt', format='%d %d %.12f\n')
##     mtx_b.save('b.txt', format='%d %d %.12f\n')
##     pause()

        output('computing resonance frequencies...')
        tt = [0]

        if isinstance(mtx_a, sc.sparse.spmatrix):
            mtx_a = mtx_a.toarray()
        if isinstance(mtx_m, sc.sparse.spmatrix):
            mtx_m = mtx_m.toarray()

        eigs, mtx_s_phi = eig(mtx_a, mtx_m, return_time=tt,
                              method=self.app_options.eigensolver)
        eigs[eigs<0.0] = 0.0
        output('...done in %.2f s' % tt[0])
        output('original eigenfrequencies:')
        output(eigs)
        opts = self.app_options
        epsilon2 = opts.scale_epsilon * opts.scale_epsilon
        eigs_rescaled = (opts.elasticity_contrast / epsilon2)  * eigs
        output('rescaled eigenfrequencies:')
        output(eigs_rescaled)
        output('number of eigenfrequencies: %d' % eigs.shape[0])

        try:
            assert_(nm.isfinite(eigs).all())
        except ValueError:
            debug()

        # B-orthogonality check.
##         print nm.dot(mtx_s_phi[:,5], nm.dot(mtx_m, mtx_s_phi[:,5]))
##         print nm.dot(mtx_s_phi[:,5], nm.dot(mtx_m, mtx_s_phi[:,0]))
##         debug()

        n_eigs = eigs.shape[0]

        variables = problem.get_variables()

        mtx_phi = nm.empty((variables.di.ptr[-1], mtx_s_phi.shape[1]),
                           dtype=nm.float64)

        make_full = variables.make_full_vec
        if eig_problem in ['simple', 'simple_liquid']:
            for ii in xrange(n_eigs):
                mtx_phi[:,ii] = make_full(mtx_s_phi[:,ii])
            eig_vectors = mtx_phi

        elif eig_problem == 'schur':
            # Update also eliminated variables.
            schur = self.app_options.schur
            primary_var = schur['primary_var']
            eliminated_var = schur['eliminated_var']

            mtx_s_phi_schur = - sc.dot(mtx_dib, mtx_s_phi)
            aux = nm.empty((variables.adi.ptr[-1],), dtype=nm.float64)
            set = variables.set_state_part
            for ii in xrange(n_eigs):
                set(aux, mtx_s_phi[:,ii], primary_var, stripped=True)
                set(aux, mtx_s_phi_schur[:,ii], eliminated_var,
                    stripped=True)

                mtx_phi[:,ii] = make_full(aux)

            indx = variables.get_indx(primary_var)
            eig_vectors = mtx_phi[indx,:]

        save = self.app_options.save
        out = {}
        state = problem.create_state()
        for ii in xrange(n_eigs):
            if (ii >= save[0]) and (ii < (n_eigs - save[1])): continue
            state.set_full(mtx_phi[:,ii], force=True)
            aux = state.create_output_dict()
            for name, val in aux.iteritems():
                out[name+'%03d' % ii] = val

        if post_process_hook is not None:
            out = post_process_hook(out, problem, mtx_phi)

        problem.domain.mesh.write(ofn_trunk + '.vtk', io='auto', out=out)

        fd = open(ofn_trunk + '_eigs.txt', 'w')
        eigs.tofile(fd, ' ')
        fd.close()

        evp = Struct(kind=eig_problem,
                     eigs=eigs, eigs_rescaled=eigs_rescaled,
                     eig_vectors=eig_vectors)
        self.cached_evp = evp

        return evp

    def eval_homogenized_coefs(self):
        if self.cached_coefs is not None:
            return self.cached_coefs

        opts = self.app_options

        if opts.homogeneous:
            rtm = opts.region_to_material
            mat_region = rtm.keys()[0]
            mat_name = rtm[mat_region]

            self.problem.update_materials()

            mat = self.problem.materials[mat_name]
            coefs = mat.get_data(mat_region, 0, opts.tensor_names)

        else:
            dc = opts.dispersion_conf
            dconf = ProblemConf.from_dict(dc['input'], dc['module'])

            dconf.materials = self.conf.materials
            dconf.regions.update(self.conf.regions)
            dconf.options['output_dir'] = self.problem.output_dir

            volume = opts.volume(self.problem, 'Y')
            problem = ProblemDefinition.from_conf(dconf, init_equations=False)
            he = HomogenizationEngine(problem, self.options, volume=volume)
            coefs = he()

##         print coefs
##         pause()
        output.prefix = self.output_prefix

        self.cached_coefs = coefs

        return coefs

    def compute_cat(self, ret_iw_dir=False):
        """Compute the Christoffel acoustic tensor, given the incident wave
        direction."""
        opts = self.app_options
        iw_dir = nm.array(opts.incident_wave_dir, dtype=nm.float64)

        dim = self.problem.get_dim()
        assert_(dim == iw_dir.shape[0])

        iw_dir = iw_dir / nla.norm(iw_dir)

        if self.cached_christoffel is not None:
            christoffel = self.cached_christoffel

        else:
            coefs = self.eval_homogenized_coefs()
            christoffel = compute_cat(coefs, iw_dir,
                                      self.app_options.dispersion)
            report_iw_cat(iw_dir, christoffel)

            self.cached_christoffel = christoffel

        if ret_iw_dir:
            return christoffel, iw_dir

        else:
            return christoffel

    def compute_phase_velocity(self):
        from sfepy.homogenization.phono import compute_density_volume_info
        opts = self.app_options
        dim = self.problem.domain.mesh.dim

        christoffel = self.compute_cat()

        self.problem.update_materials()
        dv_info = compute_density_volume_info(self.problem, opts.volume,
                                              opts.region_to_material)
        output('average density:', dv_info.average_density)

        eye = nm.eye(dim, dim, dtype = nm.float64)
        mtx_mass = eye * dv_info.average_density

        meigs, mvecs = eig(mtx_mass, mtx_b=christoffel,
                           eigenvectors=True, method=opts.eigensolver)
        phase_velocity = 1.0 / nm.sqrt(meigs)

        return phase_velocity
