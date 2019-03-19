"""
Multiprocessing functions.
"""
try:
    from mpi4py import MPI
    use_multiprocessing_mpi = MPI.COMM_WORLD.Get_size() > 1
except ImportError:
    use_multiprocessing_mpi = False

try:
    import sys
    if sys.platform.startswith('win'):
        #
        # Multiprocessing-proc implementation is currently broken on Windows platform
        # TBD: Fix (engine.py ??)
        #
        use_multiprocessing_proc = False
    else:
        from multiprocessing import cpu_count
        use_multiprocessing_proc = cpu_count() > 1
except:
    use_multiprocessing_proc = False

if use_multiprocessing_mpi:
    import sfepy.base.multiproc_mpi as multiproc_mpi
if use_multiprocessing_proc:
    import sfepy.base.multiproc_proc as multiproc_proc

use_multiprocessing = use_multiprocessing_mpi or use_multiprocessing_proc

multiprocessing_mode = None
multiprocessing_module = None


def get_multiproc(mpi=False):
    global multiprocessing_mode, multiprocessing_module
    try_proc = True
    multiproc, mode = None, None
    if mpi and use_multiprocessing_mpi:
        multiproc, mode = multiproc_mpi, 'mpi'
        try_proc = False

    if try_proc and use_multiprocessing_proc:
        multiproc, mode = multiproc_proc, 'proc'

    multiprocessing_mode = mode
    multiprocessing_module = multiproc

    return multiproc, mode


def get_num_workers():
    """Get the number of slave nodes."""
    mpi = multiprocessing_mode == 'mpi'
    return multiproc_mpi.cpu_count() - 1 if mpi else\
        multiproc_proc.cpu_count()


def is_remote_dict(d):
    if multiprocessing_module is not None:
        return multiprocessing_module.is_remote_dict(d)
    else:
        return False
