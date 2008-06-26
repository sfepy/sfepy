import threading
import atexit

##
# Extended from http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/65222.
# 22.03.2007, c
class TaskThread( threading.Thread ):
    """Thread that executes a given function every N seconds"""
    
    def __init__( self, fun, args = (), interval = 1.0 ):
        """
        Create the task thread with:
        fun      .. function to run
        args     .. tuple of function arguments
        interval .. number of seconds we sleep between executing the function
        """   
        threading.Thread.__init__( self )

        self.finished = threading.Event()
        self.interval = interval
        self.fun = fun
        self.args = args
        self.canRun = False

        # Respond to Ctrl_C.
        atexit.register( self.shutdown )
        
    def setInterval( self, interval ):
        """Set the number of seconds we sleep between executing the function."""
        self.interval = interval
    
    def shutdown( self ):
        """Stop this thread"""
        self.finished.set()
    
    def run(self):
        while 1:
            if self.finished.isSet(): return

            if self.canRun:
                self.task()
            
            # sleep for interval or until shutdown
            self.finished.wait( self.interval )
    
    def task(self):
        apply( self.fun, self.args )

    def on( self ):
        self.canRun = True

    def off( self ):
        self.canRun = False

if __name__ == '__main__':
    import pylab

    # Make a figure interactively, so that a window appears.
    pylab.ion()
    fig = pylab.figure()
    pylab.ioff()

    # Start the redrawing thread. The redrawing is switched off.
    thread = TaskThread( fig.canvas.draw )
    thread.start()

    # Switch off the redrawing (no-operation here...)
    thread.off()
    # Make some plots.
    ax = fig.add_subplot( 111 )
    ax.plot( range( 10 ) )
    # Switch on the redrawing.
    thread.on()

    # Do some work.
    for ii in xrange( 10000000 ):
        if not( ii % 1e5 ):
            print 'heavy work!', ii

    # Finish the thread.
    thread.shutdown()
    thread.join()

    # Show() the figure.
    pylab.show()
