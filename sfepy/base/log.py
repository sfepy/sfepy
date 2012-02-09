import time
import os
import atexit

import numpy as nm

from sfepy.base.base import sfepy_config_dir, ordered_iteritems
from sfepy.base.base import output, get_default, set_defaults, Output, Struct
from sfepy.base.tasks import Process, Pipe
from sfepy.linalg import cycle

_msg_no_live = """warning: log plot is disabled, install matplotlib
         (use GTKAgg backend) and multiprocessing"""

def get_logging_conf(conf):
    """
    Check for a logging configuration ('log' attribute) in conf. Supply default
    values if necessary.
    
    Returns
    -------
    log : dict
        The dictionary {'plot' : <figure_file>, 'text' : <text_log_file>}. One
        or both values can be None.
    """
    get = conf.get_default_attr

    log = get('log', None)

    default_log = {'text' : None, 'plot' : None}

    if log is None:
        log = default_log

    else:
        set_defaults(log, default_log)

    return log

class ProcessPlotter( Struct ):
    output = Output('plotter:')
    output = staticmethod( output )

    def __init__( self, aggregate = 100 ):
        Struct.__init__( self, aggregate = aggregate )

    def process_command( self, command ):
        from matplotlib.ticker import LogLocator, AutoLocator

        self.output( command[0] )

        if command[0] == 'ig':
            self.ig = command[1]

        elif command[0] == 'plot':
            xdata, ydata = command[1:]

            ig = self.ig
            ax = self.ax[ig]
            ax.set_yscale(self.yscales[ig])
            ax.yaxis.grid(True)
            ax.plot(xdata, ydata)

            if self.yscales[ig] == 'log':
                ymajor_formatter = ax.yaxis.get_major_formatter()
                ymajor_formatter.label_minor(True)
                yminor_locator = LogLocator()
            else:
                yminor_locator = AutoLocator()
                self.ax[ig].yaxis.set_minor_locator(yminor_locator)

        elif command[0] == 'vline':
            x, kwargs = command[1:]

            self.vlines[self.ig].append((x, kwargs))

        elif command[0] == 'clear':
            self.ax[self.ig].cla()

        elif command[0] == 'legends':
            for ig, ax in enumerate(self.ax):
                try:
                    ax.legend(self.data_names[ig])
                except:
                    pass

                if self.xlabels[ig]:
                    ax.set_xlabel(self.xlabels[ig])
                if self.ylabels[ig]:
                    ax.set_ylabel(self.ylabels[ig])

                for x, kwargs in self.vlines[ig]:
                    ax.axvline(x, **kwargs)

        elif command[0] == 'add_axis':
            ig, names, yscale, xlabel, ylabel = command[1:]
            self.data_names[ig] = names
            self.yscales[ig] = yscale
            self.xlabels[ig] = xlabel
            self.ylabels[ig] = ylabel
            self.n_gr = len(self.data_names)
            self.make_axes()

        elif command[0] == 'save':
            self.fig.savefig(command[1])
            self.pipe.send(True) # Acknowledge save.

    def terminate( self ):
        if self.ii:
            self.output( 'processed %d commands' % self.ii )
        self.output( 'ended.' )
        self.plt.close( 'all' )

    def poll_draw( self ):

        def call_back():
            self.ii = 0
            
            while 1:
                if not self.pipe.poll():
                    break

                command = self.pipe.recv()
                can_break = False

                if command is None:
                    self.terminate()
                    return False
                elif command[0] == 'continue':
                    can_break = True
                else:
                    self.process_command( command )

                if (self.ii >= self.aggregate) and can_break:
                    break

                self.ii += 1

            if self.ii:
                self.fig.canvas.draw()
                self.output( 'processed %d commands' % self.ii )

            return True

        return call_back

    def make_axes(self):
        self.fig.clf()
        self.ax = []

        n_col = min(5.0, nm.fix(nm.sqrt(self.n_gr)))
        if int(n_col) == 0:
            n_row = 0
        else:
            n_row = int(nm.ceil(self.n_gr / n_col))
            n_col = int(n_col)

        for ii, (ir, ic) in enumerate(cycle((n_col, n_row))):
            if ii == self.n_gr: break
            self.ax.append(self.fig.add_subplot(n_row, n_col, ii+1))
            self.vlines.setdefault(ii, [])
    
    def __call__(self, pipe, log_file, data_names, yscales, xlabels, ylabels):
        """Sets-up the plotting window, sets GTK event loop timer callback to
        callback() returned by self.poll_draw(). The callback does the actual
        plotting, taking commands out of `pipe`, and is called every second.

        Note that pyplot _must_ be imported here and not in this module so that
        the import occurs _after_ the plotting process is started in that
        process.
        """
        import gobject
        import matplotlib.pyplot as plt
        self.plt = plt

        self.output.set_output(filename=log_file)
        self.output( 'starting plotter...' )

        self.pipe = pipe
        self.data_names = data_names
        self.yscales = yscales
        self.xlabels = xlabels
        self.ylabels = ylabels
        self.n_gr = len(data_names)
        self.vlines = {}

        self.fig = self.plt.figure()
        self.make_axes()
        self.gid = gobject.timeout_add( 1000, self.poll_draw() )

        self.output( '...done' )
        self.plt.show()

def name_to_key( name, ii ):
    return name + (':%d' % ii)

class Log( Struct ):
    """Log data and (optionally) plot them in the second process via
    ProcessPlotter."""
    count = -1

    def from_conf( conf, data_names ):
        """`data_names` ... tuple of names grouped by subplots:
                            ([name1, name2, ...], [name3, name4, ...], ...)
        where name<n> are strings to display in (sub)plot legends."""
        if not isinstance( data_names, tuple ):
            data_names = (data_names,)

        obj = Log(data_names, **conf)

        return obj
    from_conf = staticmethod( from_conf )

    def __init__(self, data_names=None, yscales=None,
                 xlabels=None, ylabels=None, is_plot=True, aggregate=200,
                 formats=None, log_filename=None):
        """`data_names` ... tuple of names grouped by subplots:
                            ([name1, name2, ...], [name3, name4, ...], ...)
        where name<n> are strings to display in (sub)plot legends."""
        try:
            import matplotlib as mpl
        except:
            mpl = None

        if (mpl is not None) and mpl.rcParams['backend'] == 'GTKAgg':
            can_live_plot = True
        else:
            can_live_plot = False

        Struct.__init__(self, data_names = {},
                        n_arg = 0, n_gr = 0,
                        data = {}, x_values = {}, n_calls = 0,
                        yscales = {}, xlabels = {}, ylabels = {},
                        plot_pipe = None, formats = {}, output = None)

        if data_names is not None:
            n_gr = len(data_names)
        else:
            n_gr = 0
            data_names = []

        yscales = get_default(yscales, ['linear'] * n_gr)
        xlabels = get_default(xlabels, ['iteration'] * n_gr)
        ylabels = get_default(ylabels, [''] * n_gr )

        if formats is None:
            formats = [None] * n_gr

        for ig, names in enumerate(data_names):
            self.add_group(names, yscales[ig], xlabels[ig], ylabels[ig],
                           formats[ig])

        self.is_plot = get_default( is_plot, True )
        self.aggregate = get_default( aggregate, 100 )

        self.can_plot = (can_live_plot and (mpl is not None)
                         and (Process is not None))

        if log_filename is not None:
            self.output = Output('', filename=log_filename)
            self.output('# started: %s' % time.asctime())

        if self.is_plot and (not self.can_plot):
            output(_msg_no_live)

    def add_group(self, names, yscale=None, xlabel=None, ylabel=None,
                  formats=None):
        """Add a new data group. Notify the plotting process if it is
        already running."""
        ig = self.n_gr
        self.n_gr += 1

        self.x_values[ig] = []

        self.data_names[ig] = names
        self.yscales[ig] = yscale
        self.xlabels[ig] = xlabel
        self.ylabels[ig] = ylabel
        
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


    def __call__( self, *args, **kwargs ):
        """Log the data passed via *args, and send them to the plotting
        process, if available."""
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
            self.plot_pipe.send( ['save', save_figure] )
            ok = self.plot_pipe.recv()

        if finished:
            self.terminate()
            return

        ls = len( args ), self.n_arg
        if full and (ls[0] != ls[1]):
            if kwargs:
                return
            else:
                msg = 'log called with wrong number of arguments! (%d == %d)' \
                      % ls
                raise IndexError( msg )

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
            if isinstance( aux, nm.ndarray ):
                aux = nm.array( aux, ndmin = 1 )
                if len( aux ) == 1:
                    aux = aux[0]
                else:
                    raise ValueError, 'can log only scalars (%s)' % aux
            key = name_to_key( name, ii )
            self.data[key].append( aux )

            if self.output:
                self.output(('%%s: %%s: %s' % self.formats[key])
                            % (name, self.x_values[ig][-1], aux))

        if self.is_plot and self.can_plot:
            if self.n_calls == 0:
                atexit.register( self.terminate )

                self.__class__.count += 1

                self.plot_pipe, plotter_pipe = Pipe()
                self.plotter = ProcessPlotter( self.aggregate )
                self.plot_process = Process( target = self.plotter,
                                             args = (plotter_pipe,
                                                     self.get_log_name(),
                                                     self.data_names,
                                                     self.yscales,
                                                     self.xlabels,
                                                     self.ylabels) )
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
                for name in names:
                    key = name_to_key(name, ii)
                    try:
                        send(['plot',
                              nm.array(self.x_values[ig]),
                              nm.array(self.data[key])])
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
        """Plot vertical lines in axes given by igs at current x locations
        to mark some events."""
        if self.plot_pipe is not None:
            send = self.plot_pipe.send

            if igs is None:
                igs = range(self.n_gr)

            for ig in igs:
                x = self.x_values[ig]
                if len(x):
                    send(['ig', ig])
                    send(['vline', x[-1], kwargs])

            send(['continue'])

        if self.output:
            self.output('-----')
