"""
Functions for setting up interactive sessions.
"""
from sfepy.base.base import output, sfepy_config_dir

preexec_source = """\
from sfepy.base.base import *
from sfepy.fem import *
from sfepy.applications import pde_solve
try:
    import matplotlib as mpl
    mpl.use('WXAgg')
    from matplotlib.pyplot import *
except:
    pass
"""

preexec_source_viewer = """\
from sfepy.postprocess import Viewer
"""

no_ipython = """\
Couldn't locate IPython. Having IPython installed is greatly recommended. See
http://ipython.scipy.org for more details. You can try to install the 'ipython'
package and start isfepy again.
"""

startup_message = """
These commands were executed:

%s
Basic Usage
-----------

The function `pde_solve` can be used to run examples in problem
description files.

When in SfePy source directory, try:
>>> pb, state = pde_solve('examples/diffusion/poisson.py')
>>> print state.get_parts()
>>> view = Viewer(pb.get_output_name())
>>> view()

When in another directory (and SfePy is installed), try:
>>> from sfepy import data_dir
>>> pb, state = pde_solve(data_dir + '/examples/diffusion/poisson.py')
>>> print state.get_parts()
>>> view = Viewer(pb.get_output_name())
>>> view()

Note: try IPython's doctest mode when pasting the above code - type
%%doctest_mode in the shell.

Advanced Usage
--------------

For more advanced use refer to "Interactive Example: Linear Elasticity"
section of the tutorial, see [1].

Check also the SfePy web site [2].

[1] http://docs.sfepy.org/doc/
[2] http://sfepy.org
"""

def _make_message(ipython=True, quiet=False, source=None):
    """
    Create a banner for an interactive session.
    """
    import sys
    from sfepy import __version__ as sfepy_version

    python_version = '%d.%d.%d' % sys.version_info[:3]

    if ipython:
        import IPython
        shell_name = 'IPython'
        python_version +=  ', Ipython %s' % IPython.__version__

    else:
        shell_name = 'Python'

    args = shell_name, sfepy_version, python_version
    message = '%s console for SfePy %s (Python %s)\n' % args

    if not quiet:
        if source is None:
            source = preexec_source

        _source = ''

        for line in source.split('\n')[:-1]:
            if not line:
                _source += '\n'

            else:
                _source += '>>> ' + line + '\n'

        message += startup_message % _source

    return message

def _init_ipython_session(is_wx, is_viewer, argv=[]):
    """
    Construct new IPython session.
    """
    import IPython
    if IPython.__version__ >= '0.11':
        from IPython.frontend.terminal import ipapp
        # use an app to parse the command line, and init config
        app = ipapp.TerminalIPythonApp()
        # don't draw IPython banner during initialization:
        app.display_banner = False
        app.initialize(argv)
        ip = app.shell

        if is_wx:
            import wx
            from IPython.lib.inputhook import enable_wx
            wxapp = wx.GetApp()
            enable_wx(wxapp)

    else:
        if is_wx:
            from IPython.Shell import IPShellWX
            ip = IPShellWX(argv)

        else:
            ip = IPython.Shell.make_IPython(argv)

    return ip

def _init_python_session():
    """
    Construct new Python session.
    """
    from code import InteractiveConsole

    class SfePyConsole(InteractiveConsole):
        """
        An interactive console with readline support.
        """

        def __init__(self):
            InteractiveConsole.__init__(self)

            try:
                import readline

            except ImportError:
                pass

            else:
                import os, atexit

                readline.parse_and_bind('tab: complete')

                if hasattr(readline, 'read_history_file'):
                    history = os.path.join(sfepy_config_dir, 'history')

                    try:
                        readline.read_history_file(history)

                    except IOError:
                        pass

                    atexit.register(readline.write_history_file, history)

    return SfePyConsole()

def init_session(ipython=None, message=None, quiet=False, silent=False,
                 is_viewer=True, is_wx=True, argv=[]):
    """
    Initialize embedded IPython or Python session.
    """
    import os, sys, atexit

    in_ipython = False

    if ipython is False:
        ip = _init_python_session()
        mainloop = ip.interact

    else:
        try:
            import IPython

        except ImportError:
            if ipython is not True:
                print no_ipython
                ip = _init_python_session()
                mainloop = ip.interact

            else:
                raise RuntimeError('IPython is not available on this system')

        else:
            ipython = True
            if IPython.__version__ >= '0.11':
                try:
                    ip = get_ipython()
                except NameError:
                    ip = None

            else:
                ip = IPython.ipapi.get()
                if ip:
                    ip = ip.IP

            if ip is not None:
                in_ipython = True

            else:
                ip = _init_ipython_session(is_wx, is_viewer, argv)

            if IPython.__version__ >= '0.11':
                # runsource is gone, use run_cell instead, which doesn't
                # take a symbol arg.  The second arg is `store_history`,
                # and False means don't add the line to IPython's history.
                ip.runsource = lambda src, symbol='exec': \
                               ip.run_cell(src, False)
                mainloop = ip.mainloop

            else:
                mainloop = ip.interact

    _preexec_source = preexec_source
    if is_viewer:
        _preexec_source += preexec_source_viewer

    ip.runsource(_preexec_source, symbol='exec')

    message = _make_message(ipython, quiet, _preexec_source)

    output.set_output(filename=os.path.join(sfepy_config_dir, 'isfepy.log'),
                      combined=silent == False,
                      append=True)
    atexit.register(output, 'isfepy finished\n' + '*' * 55)

    if not in_ipython:
        mainloop(message)
        sys.exit('Exiting ...')

    else:
        ip.write(message)
        ip.set_hook('shutdown_hook', lambda ip: ip.write('Exiting ...\n'))
