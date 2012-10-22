"""
Plotting class to be used by Log.
"""
import numpy as nm

from sfepy.base.base import Output, Struct

class LogPlotter(Struct):
    """
    LogPlotter to be used by :class:`sfepy.base.log.Log`.
    """
    output = Output('plotter:')
    output = staticmethod(output)

    def __init__(self, aggregate=100):
        Struct.__init__(self, aggregate=aggregate)

    def process_command(self, command):
        from matplotlib.ticker import LogLocator, AutoLocator

        self.output(command[0])

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

    def terminate(self):
        if self.ii:
            self.output('processed %d commands' % self.ii)
        self.output('ended.')
        self.plt.close('all')

    def poll_draw(self):

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
                    self.process_command(command)

                if (self.ii >= self.aggregate) and can_break:
                    break

                self.ii += 1

            if self.ii:
                self.fig.canvas.draw()
                self.output('processed %d commands' % self.ii)

            return True

        return call_back

    def make_axes(self):
        from sfepy.linalg import cycle

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
            self.ax.append(self.fig.add_subplot(n_row, n_col, ii + 1))
            self.vlines.setdefault(ii, [])

    def __call__(self, pipe, log_file, data_names, yscales, xlabels, ylabels):
        """
        Sets-up the plotting window, sets GTK event loop timer callback to
        callback() returned by self.poll_draw(). The callback does the actual
        plotting, taking commands out of `pipe`, and is called every second.

        Note that pyplot _must_ be imported here and not in this module so that
        the import occurs _after_ the plotting process is started in that
        process.
        """
        import matplotlib.pyplot as plt
        import gobject
        self.plt = plt

        self.output.set_output(filename=log_file)
        self.output('starting plotter...')

        self.pipe = pipe
        self.data_names = data_names
        self.yscales = yscales
        self.xlabels = xlabels
        self.ylabels = ylabels
        self.n_gr = len(data_names)
        self.vlines = {}

        self.fig = self.plt.figure()
        self.make_axes()
        self.gid = gobject.timeout_add(1000, self.poll_draw())

        self.output('...done')
        self.plt.show()
