"""
Multiprocessing functions.
"""
import sys

try:
    #
    # Multiprocessing_proc implementation is currently broken on platforms
    # using 'spawn' method as default.
    #
    # ToDo: we really need a real fix (Linux fork() method seems to be
    #       insecure/deprecated) !
    #
    if sys.platform.startswith('win'):
        use_multiprocessing = False
    elif sys.platform == 'darwin' and \
        sys.version_info.major == 3 and sys.version_info.minor == 8:
        use_multiprocessing = False
    else:
        from multiprocessing import (Queue, Manager, cpu_count,
            current_process, set_start_method)
        from concurrent.futures import ProcessPoolExecutor

        if not current_process().name.startswith('HomogWorkerProcess'):
        #     # set_start_method('forkserver')
        #     # set_forkserver_preload(["numpy", "sfepy"])
        #     # set_start_method('fork')
            set_start_method('spawn')
            use_multiprocessing = cpu_count() > 1
        else:
           use_multiprocessing = False
except:
    use_multiprocessing = False

active_workers = {}
multiproc_manager = Manager()
multiproc_dict = multiproc_manager.dict()

if use_multiprocessing:
    active_workers.update({
        'pool': None,
    })


def get_dict(name, clear=False):
    if name in multiproc_dict:
        out = multiproc_dict[name]
        if clear:
            out.clear()
    else:
        out = multiproc_dict[name] = multiproc_manager.dict()

    return out


def get_workers():
    return active_workers['pool']


def get_proc_id():
    return current_process().name


def init_workers(wdir, max_workers=None, **kwargs):
    if use_multiprocessing:
        sys.path.append(wdir)
        active_workers['pool'] = ProcessPoolExecutor(max_workers, **kwargs)

        return active_workers['pool']
