from base import *
from sfepy.base.tasks import Process, Queue
from sfepy.base.plotutils import pylab

if pylab:
    import gobject
    from matplotlib.ticker import LogLocator, AutoLocator, LogFormatter
else:
    LogFormatter = object

class MyLogFormatter( LogFormatter ):
    def __call__( self, x, pos = None ):
        self.verify_intervals()
        d = abs(self.view_interval.span())
        b=self._base
        # only label the decades
        fx = nm.log(x)/nm.log(b)
        is_decade = self.is_decade(fx)
        if not is_decade and self.label_only_base: s = ''
        elif x>10000: s= '%1.0e'%x
        elif x<1: s =  '%1.0e'%x
        else: s =  self.pprint_val(x,d)

        print d, fx, is_decade, s
        
        pause()
        return s

class ProcessPlotter( Struct ):

    def poll_draw( self ):

        def call_back():
#            gobject.timeout_add( 0, None )

            command = self.queue.get()
            print command[0]

            if command is None:
                pylab.close( 'all' )
                output( 'ended.' )
                return False

            elif command[0] == 'axes':
                self.ca = self.ax[command[1]]
            elif command[0] == 'plot':
                self.ca.cla()
                self.ca.plot( command[1] )
                
            self.fig.canvas.draw()
 #           gobject.timeout_add( 100, self.poll_draw() )
            return True

        return call_back
    
    def __call__( self, queue, n_gr ):
        set_output_prefix( 'plotter:' )
        output( 'starting plotter...' )

        self.queue = queue
        self.n_gr = n_gr


        self.fig = pylab.figure()

        self.ax = []
        for ig in range( self.n_gr ):
            isub = 100 * self.n_gr + 11 + ig
            self.ax.append( self.fig.add_subplot( isub ) )
        self.ca = self.ax[0]
        
        gobject.timeout_add( 100, self.poll_draw() )

        output( '...done' )

        pylab.show()

class Log( Struct ):

    def from_conf( conf, data_names ):
        obj = Log( data_names = data_names, seq_data_names = [], igs = [],
                   data = {}, yscales = [], n_calls = 0 )
        
        for ig, names in enumerate( obj.data_names ):
            for name in names:
                obj.data[name] = []
                obj.igs.append( ig )
                obj.seq_data_names.append( name )
                obj.yscales.append( 'linear' ) # Defaults.
        obj.n_arg = len( obj.igs )
        obj.n_gr = len( obj.data_names )
        
        try:
            obj.is_plot = conf.is_plot
        except:
            obj.is_plot = True

        try:
            obj.yscales = conf.yscales
        except:
            pass
        
        return obj
    from_conf = staticmethod( from_conf )
    
    def __call__( self, *args, **kwargs ):
        finished = False
        if kwargs:
            if kwargs.has_key( 'finished' ):
                finished = kwargs['finished']

        yscales = self.yscales
                
        ls = len( args ), self.n_arg
        if ls[0] != ls[1]:
            raise IndexError, '%d == %d' % ls

        for ii, name in enumerate( self.seq_data_names ):
            aux = args[ii]
            if isinstance( aux, nm.ndarray ):
                aux = nm.array( aux, ndmin = 1 )
                if len( aux ) == 1:
                    aux = aux[0]
                else:
                    raise ValueError, 'can log only scalars (%s)' % aux
            self.data[name].append( aux )

        if self.is_plot and (pylab is not None):
            if self.n_calls == 0:
                self.plot_queue = Queue()
                self.plotter = ProcessPlotter()
                self.plot_process = Process( target = self.plotter,
                                             args = (self.plot_queue,
                                                     self.n_gr) )
                self.plot_process.start()

            self.plot_data()
            
        self.n_calls += 1

        print finished
        if finished and self.is_plot and (pylab is not None):
            self.plot_queue.put( None )
            self.plot_process.join()

    def plot_data( self ):

        put =  self.plot_queue.put
        for ii, name in enumerate( self.seq_data_names ):
            try:
                put( ['axes', self.igs[ii]] )
                put( ['plot', self.data[name]] )
            except:
                put( None )
                self.plot_process.join()
                print ii, name, self.data[name]
                raise

##         def aux():
##             for ii, name in enumerate( self.seq_data_names ):
##                 try:
##                     ax = self.ax[self.igs[ii]]
##                     ax.set_yscale( yscales[ii] )
##                     ax.yaxis.grid( True )
##                     ax.plot( self.data[name] )

##                     if yscales[ii] == 'log':
##                         ymajor_formatter = ax.yaxis.get_major_formatter()
##                         ymajor_formatter.label_minor( True )
##                         yminor_locator = LogLocator()
##                     else:
##                         yminor_locator = AutoLocator()
##                     self.ax[ig].yaxis.set_minor_locator( yminor_locator )
## #                    self.ax[ig].yaxis.set_minor_formatter( yminor_formatter )
## #                    pylab.axis( 'image' )
## ##                     bb = [nm.min( self.data[name] ), nm.max( self.data[name] )]
## ##                     ylim = [0.99 * bb[0], 1.01 * bb[1]]
## ##                     print ylim
## ##                     ax.set_ylim( ylim )
## ##                     ax.yaxis.set_ticks( bb )
##                 except:
##                     print ii, name, self.data[name]
##                     raise
##                     pass
##             for ig in range( self.n_gr ):
##                 self.ax[ig].legend( self.data_names[ig] )
