"""
Multiprocessing functions.
"""
try:
    import sys
    #
    # Multiprocessing_proc implementation is currently broken on platforms
    # using 'spawn' method as default.
    #
    # ToDo: we really need a real fix (Linux fork() method seems to be
    #       insecure/deprecated) !
    #
    if sys.platform.startswith('win'):
        use_multiprocessing_proc = False
    elif sys.platform == 'darwin' and \
        sys.version_info.major == 3 and sys.version_info.minor == 8:
        use_multiprocessing_proc = False
    else:
        from multiprocessing import cpu_count
        use_multiprocessing_proc = cpu_count() > 1
except:
    use_multiprocessing_proc = False

if use_multiprocessing_proc:
    import sfepy.base.multiproc_proc as multiproc_proc

use_multiprocessing = use_multiprocessing_proc

multiprocessing_mode = None
multiprocessing_module = None


def get_multiproc(mpi=False):
    global multiprocessing_mode, multiprocessing_module
    try_proc = True
    multiproc, mode = None, None

    if try_proc and use_multiprocessing_proc:
        multiproc, mode = multiproc_proc, 'proc'

    multiprocessing_mode = mode
    multiprocessing_module = multiproc

    return multiproc, mode


def get_num_workers():
    """Get the number of slave nodes."""
    return multiproc_proc.cpu_count()


def is_remote_dict(d):
    if multiprocessing_module is not None:
        return multiprocessing_module.is_remote_dict(d)
    else:
        return False
