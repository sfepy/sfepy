"""
Multiprocessing functions.
"""
import sys

try:
    from multiprocessing import (Manager, cpu_count,
        current_process, set_start_method)
    from concurrent.futures import ProcessPoolExecutor

    if not current_process().name.startswith('HomogeniationWorkerProcess'):
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
