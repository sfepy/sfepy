import os.path as op
import shutil

import numpy as nm

from sfepy.base.base import ordered_iteritems, output, set_defaults, assert_
from sfepy.base.base import Struct
from sfepy.homogenization.engine import HomogenizationEngine
from sfepy.homogenization.homogen_app import HomogenizationApp
from sfepy.homogenization.coefficients import Coefficients
from sfepy.homogenization.coefs_base import CoefDummy
from sfepy.applications import PDESolverApp
from sfepy.base.plotutils import plt

def try_set_defaults(obj, attr, defaults, recur=False):
    try:
        values = getattr(obj, attr)

    except:
        values = defaults

    else:
        if recur and isinstance(values, dict):
            for key, val in values.items():
                set_defaults(val, defaults)

        else:
            set_defaults(values, defaults)

    return values

def save_raw_bg_logs(filename, logs):
    """
    Save raw band gaps `logs` into the `filename` file.
    """
    out = {}

    iranges = nm.cumsum([0] + [len(ii) for ii in logs.freqs])
    out['iranges'] = iranges
    for key, log in ordered_iteritems(logs.to_dict()):
        out[key] = nm.concatenate(log, axis=0)

    nm.savez(filename, **out)

def transform_plot_data(datas, plot_transform, conf):
    if plot_transform is not None:
        fun = conf.get_function(plot_transform[0])

    dmin, dmax = 1e+10, -1e+10
    tdatas = []
    for data in datas:
        tdata = data.copy()
        if plot_transform is not None:
            tdata = fun(tdata, *plot_transform[1:])
        dmin = min(dmin, nm.nanmin(tdata))
        dmax = max(dmax, nm.nanmax(tdata))
        tdatas.append(tdata)
    dmin, dmax = min(dmax - 1e-8, dmin), max(dmin + 1e-8, dmax)
    return (dmin, dmax), tdatas

def plot_eigs(ax, plot_rsc, plot_labels, valid, freq_range, plot_range):
    """
    Plot resonance/eigen-frequencies.

    `valid` must correspond to `freq_range`

    resonances : red
    masked resonances: dotted red
    """
    if plt is None: return
    assert_(len(valid) == len(freq_range))

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

def plot_logs(ax, plot_rsc, plot_labels,
              freqs, logs, valid, freq_range, plot_range,
              draw_eigs=True, show_legend=True):
    """
    Plot logs of min/middle/max eigs of a mass matrix.
    """
    if plt is None: return

    if draw_eigs:
        plot_eigs(ax, plot_rsc, plot_labels, valid, freq_range, plot_range)

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

    ax.set_xlabel(plot_labels['x_axis'])
    ax.set_ylabel(plot_labels['y_axis'])

    if show_legend:
        ax.legend()

def plot_gap(ax, ranges, kind, kind_desc, plot_range, plot_rsc):
    """
    Plot single band gap frequency ranges as rectangles.
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
        draw_rect(ax, ranges[0], plot_range, propagation)
    elif kind == 'w':
        draw_rect(ax, ranges[0], plot_range, weak)
    elif kind == 'wp':
        draw_rect(ax, ranges[0], plot_range, weak)
        draw_rect(ax, ranges[1], plot_range, propagation)
    elif kind == 's':
        draw_rect(ax, ranges[0], plot_range, strong)
    elif kind == 'sw':
        draw_rect(ax, ranges[0], plot_range, strong)
        draw_rect(ax, ranges[1], plot_range, weak)
    elif kind == 'swp':
        draw_rect(ax, ranges[0], plot_range, strong)
        draw_rect(ax, ranges[1], plot_range, weak)
        draw_rect(ax, ranges[2], plot_range, propagation)
    elif kind == 'is':
        draw_rect(ax, ranges[0], plot_range, strong)
    elif kind == 'iw':
        draw_rect(ax, ranges[0], plot_range, weak)
    else:
        raise ValueError('unknown band gap kind! (%s)' % kind)

def plot_gaps(ax, plot_rsc, gaps, kinds, gap_ranges, freq_range, plot_range):
    """
    Plot band gaps as rectangles.
    """
    if plt is None: return

    for ii in range(len(freq_range) - 1):
        f0, f1 = freq_range[[ii, ii+1]]
        gap = gaps[ii]
        ranges = gap_ranges[ii]

        if isinstance(gap, list):
            for ig, (gmin, gmax) in enumerate(gap):
                kind, kind_desc = kinds[ii][ig]
                plot_gap(ax, ranges[ig], kind, kind_desc, plot_range, plot_rsc)

                output(ii, gmin[0], gmax[0], '%.8f' % f0, '%.8f' % f1)
                output(' -> %s\n    %s' %(kind_desc, ranges[ig]))
        else:
            gmin, gmax = gap
            kind, kind_desc = kinds[ii]
            plot_gap(ax, ranges, kind, kind_desc, plot_range, plot_rsc)
            output(ii, gmin[0], gmax[0], '%.8f' % f0, '%.8f' % f1)
            output(' -> %s\n    %s' %(kind_desc, ranges))

def _get_fig_name(output_dir, fig_name, key, common, fig_suffix):
    """
    Construct the complete name of figure file.
    """
    name = key.replace(common, '')
    if name and (not name.startswith('_')):
        name = '_' + name

    fig_name = fig_name + name + fig_suffix
    return op.join(output_dir, fig_name)

class AcousticBandGapsApp(HomogenizationApp):
    """
    Application for computing acoustic band gaps.
    """

    @staticmethod
    def process_options(options):
        """
        Application options setup. Sets default values for missing
        non-compulsory options.
        """
        get = options.get

        default_plot_options = {'show' : True,'legend' : False,}

        aux = {
            'resonance' : 'eigenfrequencies',
            'masked' : 'masked eigenfrequencies',
            'eig_min' : r'min eig($M^*$)',
            'eig_mid' : r'mid eig($M^*$)',
            'eig_max' : r'max eig($M^*$)',
            'x_axis' : r'$\sqrt{\lambda}$, $\omega$',
            'y_axis' : r'eigenvalues of mass matrix $M^*$',
        }
        plot_labels = try_set_defaults(options, 'plot_labels', aux, recur=True)

        aux = {
            'resonance' : 'eigenfrequencies',
            'masked' : 'masked eigenfrequencies',
            'eig_min' : r'$\kappa$(min)',
            'eig_mid' : r'$\kappa$(mid)',
            'eig_max' : r'$\kappa$(max)',
            'x_axis' : r'$\sqrt{\lambda}$, $\omega$',
            'y_axis' : 'polarization angles',
        }
        plot_labels_angle = try_set_defaults(options, 'plot_labels_angle', aux)

        aux = {
            'resonance' : 'eigenfrequencies',
            'masked' : 'masked eigenfrequencies',
            'eig_min' : r'wave number (min)',
            'eig_mid' : r'wave number (mid)',
            'eig_max' : r'wave number (max)',
            'x_axis' : r'$\sqrt{\lambda}$, $\omega$',
            'y_axis' : 'wave numbers',
        }
        plot_labels_wave = try_set_defaults(options, 'plot_labels_wave', aux)

        plot_rsc =  {
            'resonance' : {'linewidth' : 0.5, 'color' : 'r', 'linestyle' : '-'},
            'masked' : {'linewidth' : 0.5, 'color' : 'r', 'linestyle' : ':'},
            'x_axis' : {'linewidth' : 0.5, 'color' : 'k', 'linestyle' : '--'},
            'eig_min' : {'linewidth' : 2.0, 'color' : (0.0, 0.0, 1.0),
                         'linestyle' : ':' },
            'eig_mid' : {'linewidth' : 2.0, 'color' : (0.0, 0.0, 0.8),
                         'linestyle' : '--' },
            'eig_max' : {'linewidth' : 2.0, 'color' : (0.0, 0.0, 0.6),
                         'linestyle' : '-' },
            'strong_gap' : {'linewidth' : 0, 'facecolor' : (0.2, 0.4, 0.2)},
            'weak_gap' : {'linewidth' : 0, 'facecolor' : (0.6, 0.8, 0.6)},
            'propagation' : {'linewidth' : 0, 'facecolor' : (1, 1, 1)},
            'params' : {'axes.labelsize': 'x-large',
                        'font.size': 14,
                        'legend.fontsize': 'large',
                        'legend.loc': 'best',
                        'xtick.labelsize': 'large',
                        'ytick.labelsize': 'large',
                        'text.usetex': True},
        }
        plot_rsc = try_set_defaults(options, 'plot_rsc', plot_rsc)

        return Struct(incident_wave_dir=get('incident_wave_dir', None),

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
        get = options.get

        incident_wave_dir=get('incident_wave_dir', None,
                              'missing "incident_wave_dir" in options!')

        return Struct(incident_wave_dir=incident_wave_dir)

    def __init__(self, conf, options, output_prefix, **kwargs):
        PDESolverApp.__init__(self, conf, options, output_prefix,
                              init_equations=False)

        self.setup_options()

        if conf._filename:
            output_dir = self.problem.output_dir
            shutil.copyfile(conf._filename,
                            op.join(output_dir, op.basename(conf._filename)))

    def setup_options(self):
        HomogenizationApp.setup_options(self)

        if self.options.phase_velocity and not self.options.plot:
            process_options = AcousticBandGapsApp.process_options_pv
        else:
            process_options = AcousticBandGapsApp.process_options
        self.app_options += process_options(self.conf.options)

    def call(self):
        """
        Construct and call the homogenization engine according to options.
        """
        options = self.options

        opts = self.app_options
        conf = self.problem.conf
        coefs_name = opts.coefs
        coef_info = conf.get(opts.coefs, None,
                             'missing "%s" in problem description!'
                             % opts.coefs)

        keys = []
        if options.detect_band_gaps:
            # Compute band gaps coefficients and data.
            keys += [key for key in coef_info if key.startswith('band_gaps')]

        if options.analyze_dispersion or options.phase_velocity:

            # Insert incident wave direction to coefficients that need it.
            for key, val in coef_info.items():
                coef_opts = val.get('options', None)
                if coef_opts is None: continue

                if (('incident_wave_dir' in coef_opts)
                    and (coef_opts['incident_wave_dir'] is None)):
                    coef_opts['incident_wave_dir'] = opts.incident_wave_dir

            if options.analyze_dispersion:
                # Compute dispersion coefficients and data.
                keys += [key for key in coef_info
                        if key.startswith('dispersion')
                        or key.startswith('polarization_angles')]

            if options.phase_velocity:
                # Compute phase velocity and its requirements.
                keys += [key for key in coef_info
                        if key.startswith('phase_velocity')]

        if not keys:
            # Compute only the eigenvalue problems.
            names = [req for req in conf.get(opts.requirements, [''])
                     if req.startswith('evp')]
            coefs = {'dummy' : {'requires' : names,
                                'class' : CoefDummy,}}
            conf.coefs_dummy = coefs
            coefs_name = 'coefs_dummy'
            keys = ['dummy']

        he_options = Struct(coefs=coefs_name, requirements=opts.requirements,
                            compute_only=keys,
                            post_process_hook=self.post_process_hook,
                            multiprocessing=False)

        volumes = {}
        if hasattr(opts, 'volumes') and (opts.volumes is not None):
            volumes.update(opts.volumes)
        elif hasattr(opts, 'volume') and (opts.volume is not None):
            volumes['total'] = opts.volume
        else:
            volumes['total'] = 1.0

        he = HomogenizationEngine(self.problem, options,
                                  app_options=he_options,
                                  volumes=volumes)
        coefs = he()

        coefs = Coefficients(**coefs.to_dict())

        coefs_filename = op.join(opts.output_dir, opts.coefs_filename)
        coefs.to_file_txt(coefs_filename + '.txt',
                          opts.tex_names,
                          opts.float_format)

        bg_keys = [key for key in coefs.to_dict()
                   if key.startswith('band_gaps')
                   or key.startswith('dispersion')]
        for key in bg_keys:
            bg = coefs.get(key)
            log_save_name = bg.get('log_save_name', None)
            if log_save_name is not None:
                filename = op.join(self.problem.output_dir, log_save_name)
                bg.save_log(filename, opts.float_format, bg)

            raw_log_save_name = bg.get('raw_log_save_name', None)
            if raw_log_save_name is not None:
                filename = op.join(self.problem.output_dir, raw_log_save_name)
                save_raw_bg_logs(filename, bg.logs)

        if options.plot:
            if options.detect_band_gaps:
                self.plot_band_gaps(coefs)

            if options.analyze_dispersion:
                self.plot_dispersion(coefs)

            if opts.plot_options['show']:
                plt.show()

        if options.phase_velocity:
            keys = [key for key in coefs.to_dict()
                    if key.startswith('phase_velocity')]
            for key in keys:
                output('%s:' % key, coefs.get(key))

        return coefs

    def plot_band_gaps(self, coefs):
        opts = self.app_options

        bg_keys = [key for key in coefs.to_dict()
                   if key.startswith('band_gaps')]

        plot_opts =  opts.plot_options
        plot_rsc = opts.plot_rsc

        plt.rcParams.update(plot_rsc['params'])

        for key in bg_keys:
            bg = coefs.get(key)

            plot_labels =  opts.plot_labels.get(key, opts.plot_labels)

            plot_range, teigs = transform_plot_data(bg.logs.eigs,
                                                    opts.plot_transform,
                                                    self.conf)
            fig, ax = plt.subplots()
            plot_gaps(ax, plot_rsc, bg.gaps, bg.kinds, bg.gap_ranges,
                      bg.freq_range_margins, plot_range)
            plot_logs(ax, plot_rsc, plot_labels, bg.logs.freqs, teigs,
                      bg.valid[bg.eig_range],
                      bg.freq_range_initial,
                      plot_range,
                      show_legend=plot_opts['legend'])

            ax.set_xlim([bg.logs.freqs[0][0], bg.logs.freqs[-1][-1]])
            ax.set_ylim(plot_range)
            plt.tight_layout()

            if opts.fig_name is not None:
                fig_name = _get_fig_name(self.problem.output_dir, opts.fig_name,
                                         key, 'band_gaps', opts.fig_suffix)
                fig.savefig(fig_name)

    def plot_dispersion(self, coefs):
        opts = self.app_options

        bg_keys = [key for key in coefs.to_dict()
                   if key.startswith('dispersion')]

        plot_rsc = opts.plot_rsc
        plot_opts =  opts.plot_options
        plt.rcParams.update(plot_rsc['params'])

        plot_labels =  opts.plot_labels_angle

        for key in bg_keys:
            pas_key = key.replace('dispersion', 'polarization_angles')
            pas = coefs.get(pas_key)

            aux = transform_plot_data(pas,
                                      opts.plot_transform_angle,
                                      self.conf)
            plot_range, pas = aux


            bg = coefs.get(key)

            fig_1, ax_1 = plt.subplots()
            plot_gaps(ax_1, plot_rsc, bg.gaps, bg.kinds, bg.gap_ranges,
                      bg.freq_range_margins, plot_range)
            plot_logs(ax_1, plot_rsc, plot_labels, bg.logs.freqs, pas,
                      bg.valid[bg.eig_range],
                      bg.freq_range_initial,
                      plot_range,
                      show_legend=plot_opts['legend'])

            ax_1.set_xlim(
                [bg.logs.freqs[0][0], bg.logs.freqs[-1][-1]])
            ax_1.set_ylim(plot_range)
            plt.tight_layout()

            fig_name = opts.fig_name_angle
            if fig_name is not None:
                fig_name = _get_fig_name(self.problem.output_dir, fig_name,
                                         key, 'dispersion', opts.fig_suffix)
                fig_1.savefig(fig_name)

            aux = transform_plot_data(bg.logs.eigs,
                                      opts.plot_transform_wave,
                                      self.conf)
            plot_range, teigs = aux

            plot_labels =  opts.plot_labels_wave

            fig_2, ax_2 = plt.subplots()
            plot_gaps(ax_2, plot_rsc, bg.gaps, bg.kinds, bg.gap_ranges,
                      bg.freq_range_margins, plot_range)
            plot_logs(ax_2, plot_rsc, plot_labels, bg.logs.freqs, teigs,
                      bg.valid[bg.eig_range],
                      bg.freq_range_initial,
                      plot_range,
                      show_legend=plot_opts['legend'])

            ax_2.set_xlim([bg.logs.freqs[0][0], bg.logs.freqs[-1][-1]])
            ax_2.set_ylim(plot_range)
            plt.tight_layout()

            fig_name = opts.fig_name_wave
            if fig_name is not None:
                fig_name = _get_fig_name(self.problem.output_dir, fig_name,
                                         key, 'dispersion', opts.fig_suffix)
                fig_2.savefig(fig_name)

        if plot_opts['show']:
            plt.show()
