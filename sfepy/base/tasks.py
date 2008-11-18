import threading
import atexit

try:
    import multiprocessing
except ImportError:
    multiprocessing = None

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
        self.can_run = False

        # Respond to Ctrl_C.
        atexit.register( self.shutdown )
        
    def set_interval( self, interval ):
        """Set the number of seconds we sleep between executing the function."""
        self.interval = interval
    
    def shutdown( self ):
        """Stop this thread"""
        self.finished.set()
    
    def run(self):
        while 1:
            if self.finished.isSet(): return

            if self.can_run:
                self.task()
            
            # sleep for interval or until shutdown
            self.finished.wait( self.interval )
    
    def task(self):
        apply( self.fun, self.args )

    def on( self ):
        self.can_run = True

    def off( self ):
        self.can_run = False


if multiprocessing is not None:
    from Queue import Empty
    Process = multiprocessing.Process
    Queue = multiprocessing.Queue
