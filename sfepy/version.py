# SfePy version
__version__ = '2011.1'


def get_basic_info(version=__version__):
    """
    Return SfePy installation directory information. Append current git
    commit hash to `version`.
    """
    import os.path as op

    # If installed, up_dir is '.', otherwise (in (git) source directory) '..'.
    for up_dir in ['..', '.']:
        top_dir = op.normpath(op.join(op.dirname(__file__), up_dir))
        aux = op.join(top_dir, 'README')
        if op.isfile(aux):
            break
    else:
        raise RuntimeError('cannot determine SfePy top level directory!')

    # Append current git commit hash to __version__.
    master = op.join(top_dir, '.git/refs/heads/master')
    if op.isfile(master):
        fd = open(master, 'r')
        version += ' (%s)' % fd.readline().strip()
        fd.close()

    in_source_tree = up_dir == '..'

    return version, top_dir, in_source_tree

__version__, top_dir, in_source_tree = get_basic_info(__version__)
