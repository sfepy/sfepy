from base import *
from sfepy.base.tasks import TaskThread

try:
    import pylab
    from matplotlib.ticker import LogLocator, AutoLocator, LogFormatter
except ImportError:
    print 'matplotlib import failed!'
    pylab = None
    LogFormatter = object
    
class MyLogFormatter( LogFormatter ):
    def __call__( self, x, pos = None ):
        self.verify_intervals()
        d = abs(self.viewInterval.span())
        b=self._base
        # only label the decades
        fx = nm.log(x)/nm.log(b)
        isDecade = self.is_decade(fx)
        if not isDecade and self.labelOnlyBase: s = ''
        elif x>10000: s= '%1.0e'%x
        elif x<1: s =  '%1.0e'%x
        else: s =  self.pprint_val(x,d)

        print d, fx, isDecade, s
        
        pause()
        return s

##
# 20.04.2006, c
class Log( Struct ):

    ##
    # 20.04.2006, c
    # 21.04.2006
    # 21.03.2007
    def fromConf( conf, dataNames ):
        obj = Log( dataNames = dataNames, seqDataNames = [], igs = [],
                   data = {}, yscales = [], nCalls = 0 )
        
        for ig, names in enumerate( obj.dataNames ):
            for name in names:
                obj.data[name] = []
                obj.igs.append( ig )
                obj.seqDataNames.append( name )
                obj.yscales.append( 'linear' ) # Defaults.
        obj.nArg = len( obj.igs )
        obj.nGr = len( obj.dataNames )
        
        try:
            obj.isPlot = conf.isPlot
        except:
            obj.isPlot = True

        try:
            obj.yscales = conf.yscales
        except:
            pass
        
        return obj
    fromConf = staticmethod( fromConf )
    
    ##
    # 20.04.2006, c
    # 21.04.2006
    # 26.04.2006
    # 28.08.2006
    # 21.03.2007
    # 22.03.2007
    def __call__( self, *args, **kwargs ):
        finished = False
        if kwargs:
            if kwargs.has_key( 'finished' ):
                finished = kwargs['finished']

        yscales = self.yscales
                
        ls = len( args ), self.nArg
        if ls[0] != ls[1]:
            raise IndexError, '%d == %d' % ls

        for ii, name in enumerate( self.seqDataNames ):
            aux = args[ii]
            if isinstance( aux, nm.ndarray ):
                aux = nm.array( aux, ndmin = 1 )
                if len( aux ) == 1:
                    aux = aux[0]
                else:
                    raise ValueError, 'can log only scalars (%s)' % aux
            self.data[name].append( aux )

        if self.isPlot and (pylab is not None):
            if self.nCalls == 0:
                pylab.ion()
                self.fig = pylab.figure()
                pylab.ioff()
                self.ax = []
                for ig in range( self.nGr ):
                    isub = 100 * self.nGr + 11 + ig
                    self.ax.append( self.fig.add_subplot( isub ) )
                self.thread = TaskThread( self.fig.canvas.draw )
                self.thread.start()
            
            self.thread.off()
            for ig in range( self.nGr ):
                self.ax[ig].cla()

            for ii, name in enumerate( self.seqDataNames ):
                try:
                    ax = self.ax[self.igs[ii]]
                    ax.set_yscale( yscales[ii] )
                    ax.yaxis.grid( True )
                    ax.plot( self.data[name] )

                    if yscales[ii] == 'log':
                        ymajorFormatter = ax.yaxis.get_major_formatter()
                        ymajorFormatter.label_minor( True )
                        yminorLocator = LogLocator()
                    else:
                        yminorLocator = AutoLocator()
                    self.ax[ig].yaxis.set_minor_locator( yminorLocator )
#                    self.ax[ig].yaxis.set_minor_formatter( yminorFormatter )
#                    pylab.axis( 'image' )
##                     bb = [nm.min( self.data[name] ), nm.max( self.data[name] )]
##                     ylim = [0.99 * bb[0], 1.01 * bb[1]]
##                     print ylim
##                     ax.set_ylim( ylim )
##                     ax.yaxis.set_ticks( bb )
                except:
                    print ii, name, self.data[name]
                    raise
                    pass
            for ig in range( self.nGr ):
                self.ax[ig].legend( self.dataNames[ig] )
            self.thread.on()
#            pause()
            
        self.nCalls += 1

        if finished and self.isPlot and (pylab is not None):
            self.thread.shutdown()
            self.thread.join()
            pylab.show()

            return
