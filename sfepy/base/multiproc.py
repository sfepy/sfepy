"""
Multiprocessing functions.
"""
try:
    from mpi4py import MPI
    use_multiprocessing_mpi = MPI.COMM_WORLD.Get_size() > 1
except:
    use_multiprocessing_mpi = False

try:
    from multiprocessing import cpu_count
    use_multiprocessing_threads = cpu_count() > 1
except:
    use_multiprocessing_threads = False

if use_multiprocessing_mpi:
    import sfepy.base.multiproc_mpi as multiproc_mpi
if use_multiprocessing_threads:
    import sfepy.base.multiproc_threads as multiproc_threads

use_multiprocessing = use_multiprocessing_mpi or use_multiprocessing_threads

multiprocessing_mode = None
multiprocessing_module = None


def get_multiproc(mpi=False):
    global multiprocessing_mode, multiprocessing_module
    try_threads = True
    multiproc, mode = None, None
    if mpi and use_multiprocessing_mpi:
        multiproc, mode = multiproc_mpi, 'mpi'
        try_threads = False

    if try_threads and use_multiprocessing_threads:
        multiproc, mode = multiproc_threads, 'threads'

    multiprocessing_mode = mode
    multiprocessing_module = multiproc

    return multiproc, mode


def get_num_workers():
    """Get the number of slave nodes."""
    mpi = multiprocessing_mode == 'mpi'
    return multiproc_mpi.cpu_count() - 1 if mpi else\
        multiproc_threads.cpu_count()


def is_remote_dict(d):
    if multiprocessing_module is not None:
        return multiprocessing_module.is_remote_dict(d)
    else:
        return False
