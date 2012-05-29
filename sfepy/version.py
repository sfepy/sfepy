# SfePy version
__version__ = '2012.2'

# "Minimal" supported versions.
NUMPY_MIN_VERSION = '1.3'
SCIPY_MIN_VERSION = '0.7'
MATPLOTLIB_MIN_VERSION = '0.99.0'
PYPARSING_MIN_VERSION = '1.5.0'
PYTABLES_MIN_VERSION = '2.1.2'
IPYTHON_MIN_VERSION = '0.10.0'
MAYAVI_MIN_VERSION = '3.3.0'
SYMPY_MIN_VERSION = '0.6.7'
CYTHON_MIN_VERSION = '0.14.1'

def get_basic_info(version=__version__):
    """
    Return SfePy installation directory information. Append current git
    commit hash to `version`.
    """
    import os.path as op
    from sfepy import Config

    # If installed, up_dir is '.', otherwise (in (git) source directory) '..'.
    for up_dir in ['..', '.']:
        top_dir = op.normpath(op.join(op.dirname(__file__), up_dir))
        aux = op.join(top_dir, 'README')
        if op.isfile(aux):
            break
    else:
        raise RuntimeError('cannot determine SfePy top level directory!')

    config = Config()
    if not config.is_release():
        # Append current git commit hash to __version__.
        master = op.join(top_dir, '.git/refs/heads/master')
        if op.isfile(master):
            fd = open(master, 'r')
            version += '-git-%s' % fd.readline().strip()
            fd.close()

    in_source_tree = up_dir == '..'

    return version, top_dir, in_source_tree

__version__, top_dir, in_source_tree = get_basic_info(__version__)
