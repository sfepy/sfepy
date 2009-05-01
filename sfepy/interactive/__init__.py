from sfepy.base.base import *
from sfepy.fem import *
from sfepy.applications import pde_solve
from sfepy.postprocess import Viewer

sfepy_config_dir = os.path.expanduser('~/.sfepy')
if not os.path.exists(sfepy_config_dir):
    os.makedirs(sfepy_config_dir)

def init_session(session, message=None, silent=False, argv=[]):
    """Initialize embedded IPython or Python session. """
    import os, sys, atexit

    def init_IPython():
        return IPython.Shell.make_IPython(argv)

    def init_Python():
        import code

        class HistoryConsole(code.InteractiveConsole):
            def __init__(self):
                code.InteractiveConsole.__init__(self)

                history = os.path.join(sfepy_config_dir, 'history')
                try:
                    import readline

                    readline.parse_and_bind('tab: complete')

                    if hasattr(readline, 'read_history_file'):
                        try:
                            readline.read_history_file(history)
                        except IOError:
                            pass

                        atexit.register(readline.write_history_file, history)
                except ImportError:
                    pass

        return HistoryConsole()

    if session not in ['ipython', 'python']:
        raise ValueError("'%s' is not a valid session name" % session)

    in_ipyshell = False

    try:
        import IPython

        ip = IPython.ipapi.get()
        if ip is not None:
            if session == 'ipython':
                ip, in_ipyshell = ip.IP, True
            else:
                raise ValueError("Can't start Python shell from IPython")
        else:
            if session == 'ipython':
                ip = init_IPython()
            else:
                ip = init_Python()
    except ImportError:
        if session == 'ipython':
            raise
        else:
            ip = init_Python()

    ip.runcode(ip.compile("from sfepy.interactive import *"))

    output.set_output(filename=os.path.join(sfepy_config_dir, 'isfepy.log'),
                      combined=silent == False,
                      append=True)
    atexit.register(output, 'isfepy finished\n' + '*' * 55)

    if not in_ipyshell:
        import sfepy
        py_version = "%d.%d.%d" % sys.version_info[:3]

        welcome = "Python %s console for SfePy %s" \
                  % (py_version, sfepy.__version__)

        if message is not None:
            message = welcome + '\n\n' + message
        else:
            message = welcome + '\n'

        ip.interact(message)
        sys.exit('Exiting ...')
    else:
        def shutdown_hook(self):
            print "Exiting ..."

        ip.set_hook('shutdown_hook', shutdown_hook)
