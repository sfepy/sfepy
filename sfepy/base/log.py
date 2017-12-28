from __future__ import absolute_import
import time
import os
import atexit
import six
from six.moves import range

try:
    from multiprocessing import Process, Pipe

except ImportError:
    Process = None

import numpy as nm

from sfepy.base.base import sfepy_config_dir, ordered_iteritems
from sfepy.base.base import output, get_default, set_defaults, Output, Struct
from sfepy.base.log_plotter import LogPlotter

_msg_no_live = 'warning: log plot is disabled, install matplotlib and' \
    ' multiprocessing'

def get_logging_conf(conf, log_name='log'):
    """
    Check for a log configuration ('log' attribute by default) in
    `conf`. Supply default values if necessary.

    Parameters
    ----------
    conf : Struct
        The configuration object.
    log_name : str, optional
        The name of the log configuration attribute in `conf`.

    Returns
    -------
    log : dict
        The dictionary {'plot' : <figure_file>, 'text' : <text_log_file>}. One
        or both values can be None.
    """
    log = conf.get(log_name, None)

    default_log = {'text' : None, 'plot' : None}

    if log is None:
        log = default_log

    else:
        set_defaults(log, default_log)

    return log

def name_to_key(name, ii):
    return name + (':%d' % ii)

def read_log(filename):
    """
    Read data saved by :class:`Log` into a text file.

    Parameters
    ----------
    filename : str
        The name of a text log file.

    Returns
    -------
    log : dict
        The log with data names as keys and ``(xs, ys, vlines)`` as values.
    info : dict
        The log plot configuration with subplot numbers as keys.
    """
    log = {}
    info = {}

    fd = open(filename, 'r')

    last_xval = None
    for line in fd:
        if line[0] == '#':
            ls = line.split(':')
            if ls[0] == '# groups':
                n_gr = int(ls[1])
                for ig in range(n_gr):
                    next(fd)
                    line_info = next(fd)
                    xlabel, ylabel, yscales = line_info.split(',')

                    line_names = next(fd)
                    names = line_names.split(':')[1]

                    line_plot_kwargs = next(fd)

                    info[ig] = (xlabel.split(':')[1].strip().strip('"'),
                                ylabel.split(':')[1].strip().strip('"'),
                                yscales.split(':')[1].strip().strip('"'),
                                [name.strip().strip('"')
                                 for name in names.split(',')],
                                eval(line_plot_kwargs[19:].strip().strip('"')
                                     + ','))

            continue

        ls = line.split(':')

        key = ls[0]
        xs, ys, vlines = log.setdefault(key, ([], [], []))

        if (len(ls) == 2) and (last_xval is not None):
            vlines.append(last_xval)

        else:
            try:
                xval, yval = float(ls[1]), float(ls[2])

            except ValueError:
                continue

            xs.append(xval)
            ys.append(yval)

            last_xval = xval

    fd.close()

    for key, (xs, ys, vlines) in six.iteritems(log):
        log[key] = (nm.array(xs), nm.array(ys), nm.array(vlines))

    return log, info

def plot_log(axs, log, info, xticks=None, yticks=None, groups=None):
    """
    Plot log data returned by :func:`read_log()` into a specified figure.

    Parameters
    ----------
    axs : sequence of matplotlib.axes.Axes
        The list of axes for the log data plots.
    log : dict
        The log with data names as keys and ``(xs, ys, vlines)`` as values.
    info : dict
        The log plot configuration with subplot numbers as keys.
    xticks : list of arrays, optional
        The list of x-axis ticks (array or None) for each subplot.
    yticks : list of arrays, optional
        The list of y-axis ticks (array or None) for each subplot.
    groups : list, optional
        The list of data groups subplots. If not given, all groups are plotted.
    """
    import matplotlib.pyplot as plt

    if axs is None:
        fig = plt.figure()

    else:
        fig = None

    if groups is None:
        n_gr = len(info)
        groups = nm.arange(n_gr)

    else:
        n_gr = len(groups)

    n_col = min(5.0, nm.fix(nm.sqrt(n_gr)))
    if int(n_col) == 0:
        n_row = 0

    else:
        n_row = int(nm.ceil(n_gr / n_col))
        n_col = int(n_col)

    if xticks is None:
        xticks = [None] * n_gr

    if yticks is None:
        yticks = [None] * n_gr

    isub = 0
    for ii, (xlabel, ylabel, yscale, names, plot_kwargs) in six.iteritems(info):
        if ii not in groups: continue

        if axs is None:
            ax = fig.add_subplot(n_row, n_col, isub + 1)

        else:
            ax = axs[ii]

        ax.set_yscale(yscale)

        if xlabel:
            ax.set_xlabel(xlabel)

        if ylabel:
            ax.set_ylabel(ylabel)

        for ip, name in enumerate(names):
            xs, ys, vlines = log[name]
            ax.plot(xs, ys, label=name, **plot_kwargs[ip])

            for x in vlines:
                ax.axvline(x, color='k', alpha=0.3)

        if xticks[isub] is not None:
            ax.set_xticks(xticks[isub])

        else:
            ax.locator_params(axis='x', nbins=10)

        if yticks[isub] is not None:
            ax.set_yticks(yticks[isub])

        ax.legend(loc='best')
        isub += 1

    plt.tight_layout(pad=0.5)

class Log(Struct):
    """
    Log data and (optionally) plot them in the second process via
    LogPlotter.
    """
    count = -1

    @staticmethod
    def from_conf(conf, data_names):
        """
        Parameters
        ----------
        data_names : list of lists of str
            The data names grouped by subplots: [[name1, name2, ...], [name3,
            name4, ...], ...], where name<n> are strings to display in
            (sub)plot legends.
        """
        obj = Log(data_names, **conf)

        return obj

    def __init__(self, data_names=None, plot_kwargs=None,
                 xlabels=None, ylabels=None, yscales=None,
                 is_plot=True, aggregate=100, sleep=1.0,
                 log_filename=None, formats=None):
        """
        Parameters
        ----------
        data_names : list of lists of str
            The data names grouped by subplots: [[name1, name2, ...], [name3,
            name4, ...], ...], where name<n> are strings to display in
            (sub)plot legends.
        plot_kwargs : list of (lists of dicts) or dicts
            The keyword arguments dicts passed to plot(). For each group the
            item can be either a dict that is applied to all lines in the
            group, or a list of dicts for each line in the group.
        xlabels : list of str
            The x axis labels of subplots.
        ylabels : list of str
            The y axis labels of subplots.
        yscales : list of 'linear' or 'log'
            The y axis scales of subplots.
        is_plot : bool
            If True, try to use LogPlotter for plotting.
        aggregate : int
            The number of plotting commands to process before a redraw.
        sleep : float
            The number of seconds to sleep between polling draw commands.
        log_filename : str, optional
            If given, save log data into a log file.
        formats : list of lists of number format strings
            The print formats of data to be used in a log file, group in the
            same way as subplots.
        """
        try:
            import matplotlib as mpl
        except:
            mpl = None

        Struct.__init__(self,
                        is_plot=is_plot, aggregate=aggregate, sleep=sleep,
                        data_names={}, n_arg=0, n_gr=0,
                        data={}, x_values={}, n_calls=0, plot_kwargs={},
                        yscales={}, xlabels={}, ylabels={},
                        plot_pipe=None, formats={}, output=None)

        if data_names is not None:
            n_gr = len(data_names)
        else:
            n_gr = 0
            data_names = []

        plot_kwargs = get_default(plot_kwargs, [{}] * n_gr)
        yscales = get_default(yscales, ['linear'] * n_gr)
        xlabels = get_default(xlabels, ['iteration'] * n_gr)
        ylabels = get_default(ylabels, [''] * n_gr)

        if formats is None:
            formats = [None] * n_gr

        for ig, names in enumerate(data_names):
            self.add_group(names, plot_kwargs[ig],
                           yscales[ig], xlabels[ig], ylabels[ig],
                           formats[ig])

        self.can_plot = (mpl is not None) and (Process is not None)

        if log_filename is not None:
            self.output = Output('', filename=log_filename)
            self.output('# started: %s' % time.asctime())
            self.output('# groups: %d' % n_gr)
            for ig, names in enumerate(data_names):
                self.output('#   %d' % ig)
                self.output('#     xlabel: "%s", ylabel: "%s", yscales: "%s"'
                            % (xlabels[ig], ylabels[ig], yscales[ig]))
                self.output('#     names: "%s"' % ', '.join(names))
                self.output('#     plot_kwargs: "%s"'
                            % ', '.join('%s' % ii
                                        for ii in self.plot_kwargs[ig]))

        if self.is_plot and (not self.can_plot):
            output(_msg_no_live)

    def add_group(self, names, plot_kwargs=None,
                  yscale=None, xlabel=None, ylabel=None, formats=None):
        """
        Add a new data group. Notify the plotting process if it is
        already running.
        """
        ig = self.n_gr
        self.n_gr += 1

        self.x_values[ig] = []

        self.data_names[ig] = names
        self.yscales[ig] = yscale
        self.xlabels[ig] = xlabel
        self.ylabels[ig] = ylabel

        if isinstance(plot_kwargs, dict):
            self.plot_kwargs[ig] = [plot_kwargs] * len(names)

        else:
            self.plot_kwargs[ig] = plot_kwargs

        ii = self.n_arg
        for iseq, name in enumerate(names):
            key = name_to_key(name, ii)
            self.data[key] = []
            ii += 1

            if formats is not None:
                self.formats[key] = formats[iseq]
            else:
                self.formats[key] = '%.3e'

        self.n_arg = ii

        if self.plot_pipe is not None:
            send = self.plot_pipe.send
            send(['add_axis', ig, names, yscale, xlabel, ylabel])

    def iter_names(self, igs=None):
        if igs is None:
            igs = nm.arange(self.n_gr)

        ii = iseq = 0
        for ig, names in ordered_iteritems(self.data_names):
            for name in names:
                if ig in igs:
                    yield ig, ii, iseq, name
                    iseq += 1
                ii += 1

    def get_log_name(self):
        return os.path.join(sfepy_config_dir,
                            'plotter_%03d.log' % self.__class__.count)


    def __call__(self, *args, **kwargs):
        """
        Log the data passed via *args, and send them to the plotting
        process, if available.
        """
        finished = False
        save_figure = ''
        x_values = None
        igs = nm.arange(self.n_gr)
        full = True
        if kwargs:
            if 'finished' in kwargs:
                finished = kwargs['finished']
            if 'save_figure' in kwargs:
                save_figure = kwargs['save_figure']
            if 'x' in kwargs:
                x_values = kwargs['x']

            if 'igs' in kwargs:
                igs = nm.array(kwargs['igs'])
                full = False

        if save_figure and (self.plot_pipe is not None):
            self.plot_pipe.send(['save', save_figure])
            self.plot_pipe.recv()

        if finished:
            self.terminate()
            return

        ls = len(args), self.n_arg
        if full and (ls[0] != ls[1]):
            if kwargs:
                return
            else:
                msg = 'log called with wrong number of arguments! (%d == %d)' \
                      % ls
                raise IndexError(msg)

        for ig in igs:
            if (x_values is not None) and (x_values[ig] is not None):
                self.x_values[ig].append(x_values[ig])
            else:
                if len(self.x_values[ig]):
                    ii = self.x_values[ig][-1] + 1
                else:
                    ii = 0
                self.x_values[ig].append(ii)

        for ig, ii, iseq, name in self.iter_names(igs):
            aux = args[iseq]
            if isinstance(aux, nm.ndarray):
                aux = nm.array(aux, ndmin = 1)
                if len(aux) == 1:
                    aux = aux[0]
                else:
                    raise ValueError('can log only scalars (%s)' % aux)
            key = name_to_key(name, ii)
            self.data[key].append(aux)

            if self.output:
                self.output(('%%s: %%s: %s' % self.formats[key])
                            % (name, self.x_values[ig][-1], aux))

        if self.is_plot and self.can_plot:
            if self.n_calls == 0:
                atexit.register(self.terminate)

                self.__class__.count += 1

                self.plot_pipe, plotter_pipe = Pipe()
                self.plotter = LogPlotter(self.aggregate, self.sleep)
                self.plot_process = Process(target=self.plotter,
                                            args=(plotter_pipe,
                                                  self.get_log_name(),
                                                  self.data_names,
                                                  self.yscales,
                                                  self.xlabels,
                                                  self.ylabels))
                self.plot_process.daemon = True
                self.plot_process.start()

            self.plot_data(igs)

        self.n_calls += 1

    def terminate(self):
        if self.output is not None:
            self.output('# ended: %s' % time.asctime())
            self.output = None

        if self.is_plot and self.can_plot:
            self.plot_pipe.send(None)
            self.plot_process.join()
            self.n_calls = 0
            output('terminated')

    def plot_data(self, igs):
        send = self.plot_pipe.send

        ii = 0
        for ig, names in ordered_iteritems(self.data_names):
            if ig in igs:
                send(['ig', ig])
                send(['clear'])
                for ip, name in enumerate(names):
                    key = name_to_key(name, ii)
                    try:
                        send(['plot',
                              nm.array(self.x_values[ig]),
                              nm.array(self.data[key]),
                              self.plot_kwargs[ig][ip]])
                    except:
                        msg = "send failed! (%s, %s, %s)!" \
                              % (ii, name, self.data[key])
                        raise IOError(msg)
                    ii += 1

            else:
                ii += len(names)

        send(['legends'])
        send(['continue'])

    def plot_vlines(self, igs=None, **kwargs):
        """
        Plot vertical lines in axes given by igs at current x locations
        to mark some events.
        """
        if igs is None:
            igs = range(self.n_gr)

        if self.plot_pipe is not None:
            send = self.plot_pipe.send

            for ig in igs:
                x = self.x_values[ig]
                if len(x):
                    send(['ig', ig])
                    send(['vline', x[-1], kwargs])

            send(['continue'])

        if self.output:
            for ig in igs:
                for name in self.data_names[ig]:
                    self.output(name + ': -----')
