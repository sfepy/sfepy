"""
Multiprocessing functions - using multiprocessing (process based) module.
"""
try:
    from multiprocessing import cpu_count, Manager, Queue, Lock,\
        managers, Process, Pool
    Process;
    use_multiprocessing = cpu_count() > 1
except:
    use_multiprocessing = False
    managers = None

try:
    import queue
except ImportError:
    import Queue as queue

global_multiproc_dict = {}


class MyQueue(object):
    def __init__(self):
        self.queue = Queue()

    def get(self):
        try:
            name = self.queue.get(False, 0.01)  # get (or wait for) a task
            return name
        except queue.Empty:
            return None

    def put(self, value):
        self.queue.put(value)


def get_manager():
    """
    Get the multiprocessing manager. If not in the global cache,
    create a new instance.

    Returns
    -------
    manager : manager
        The multiprocessing manager.
    """
    if use_multiprocessing and 'manager' not in global_multiproc_dict:
        global_multiproc_dict['manager'] = Manager()

    return global_multiproc_dict['manager']


def get_mpdict_value(mode, key, clear=False):
    """
    Get the item from the global multiprocessing cache.

    Parameters
    ----------
    mode : str
        The type of the required object.
    key : immutable type
        The key of the required object.
    clear : bool
        If True, clear the dictionary or list (for modes 'dict' and 'list').

    Returns
    -------
    value : remote object
        The remote object.
    """
    m = get_manager()
    if key in global_multiproc_dict:
        if clear:
            if mode == 'dict':
                global_multiproc_dict[key].clear()
            elif mode == 'list':
                del global_multiproc_dict[key][:]
    else:
        if mode == 'dict':
            global_multiproc_dict[key] = m.dict()
        elif mode == 'list':
            global_multiproc_dict[key] = m.list()
        elif mode == 'ivalue':
            global_multiproc_dict[key] = m.Value('i', 0)
        elif mode == 'queue':
            global_multiproc_dict[key] = MyQueue()
        elif mode == 'lock':
            global_multiproc_dict[key] = Lock()

    return global_multiproc_dict[key]


def get_dict(name, clear=False, **kwargs):
    """Get the remote dictionary."""
    return get_mpdict_value('dict', 'd_' + name, clear=clear)


def get_list(name, clear=False):
    """Get the remote list."""
    return get_mpdict_value('list', 'l_' + name, clear=clear)


def get_int_value(name, val0=0):
    """Get the remote integer value."""
    out = get_mpdict_value('ivalue', 'iv_' + name)
    out.value = val0
    return out


def get_queue(name):
    """Get the global queue."""
    return get_mpdict_value('queue', 'queue_' + name)


def get_lock(name):
    """Get the global lock."""
    return get_mpdict_value('lock', 'lock_' + name)


def is_remote_dict(d):
    """Return True if 'd' is   instance."""
    return isinstance(d, managers.DictProxy)
