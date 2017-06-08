"""
MPI multiprocessing.
"""
import logging
import os

try:
    from mpi4py import MPI
    mpi_comm = MPI.COMM_WORLD
    mpi_rank = mpi_comm.Get_rank()
    mpi_status = MPI.Status()

    use_multiprocessing = mpi_comm.Get_size() >= 2
except:
    use_multiprocessing = False

global_multiproc_dict = {}
mpi_master = 0


class MPILogFile(MPI.File):
    def write(self, *args, **kwargs):
        self.Write_shared(*args, **kwargs)


class MPIFileHandler(logging.FileHandler):
    "MPI file class for logging process communication."
    def __init__(self, filename, mode=MPI.MODE_WRONLY,
                 encoding=None, delay=0, comm=MPI.COMM_WORLD):
        self.baseFilename = os.path.abspath(filename)
        self.mode = mode
        self.encoding = encoding
        self.comm = comm
        if delay:
            logging.Handler.__init__(self)
            self.stream = None
        else:
            logging.StreamHandler.__init__(self, self._open())

    def _open(self):
        stream = MPILogFile.Open(self.comm, self.baseFilename, self.mode)
        stream.Set_atomicity(True)
        return stream

    def close(self):
        if self.stream:
            self.stream.Sync()
            self.stream.Close()
            self.stream = None


def set_logging_level(log_level='info'):
    if log_level == 'debug':
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)


def get_logger(log_filename='multiproc_mpi.log'):
    """Get the MPI logger which log information into a shared file."""
    open(log_filename, 'w').close()  # empy log file

    log_id = 'master' if mpi_rank == 0 else 'slave%d' % mpi_comm.rank
    logger = logging.getLogger(log_id)
    logger.setLevel(logging.INFO)

    mh = MPIFileHandler(log_filename)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s'
                                  + '- %(name)s: %(message)s')
    mh.setFormatter(formatter)
    logger.addHandler(mh)

    return logger


logger = get_logger()


def enum(*sequential):
    enums = dict(zip(sequential, range(len(sequential))))
    reverse = dict((value, key) for key, value in enums.iteritems())
    enums['name'] = reverse
    return type('Enum', (), enums)


tags = enum('READY', 'DONE', 'START', 'CONTINUE',
            'LOCK', 'UNLOCK', 'LOCK_STATUS',
            'SET_DICT', 'SET_DICT_IMMUTABLE', 'GET_DICT', 'DICT_VAL',
            'SET_DICT_STATUS', 'GET_DICT_KEYS', 'DICT_KEYS',
            'GET_DICT_LEN', 'DICT_LEN', 'GET_DICT_IN', 'DICT_IN')


def cpu_count():
    """Get the number of MPI nodes."""
    return mpi_comm.Get_size()


def get_int_value(name, init_value=0):
    """Get the remote integer value."""
    key = '__int_%s' % name
    rdict = get_dict(key, mutable=True)
    if mpi_rank == mpi_master:
        val = RemoteInt(rdict, init_value)
    else:
        val = RemoteInt(rdict)

    return val


def get_queue(name):
    """Get the queue."""
    key = '__queue_%s' % name
    rdict = get_dict(key, mutable=True)
    if mpi_rank == mpi_master:
        val = RemoteQueue(rdict, [])
    else:
        val = RemoteQueue(rdict)

    return val


def get_dict(name, mutable=False, clear=False):
    """Get the remote dictionary."""
    if mpi_rank == mpi_master:
        if name in global_multiproc_dict:
            if clear:
                global_multiproc_dict[name].clear()
            return global_multiproc_dict[name]
        else:
            rdict = RemoteDictMaster()
            rdict.set_name(name)
            global_multiproc_dict[name] = rdict
            return rdict
    else:
        return RemoteDict(name, mutable=mutable)


class RemoteInt(object):
    """Remote intiger class, data saved in RemoteDict."""
    class IntDesc(object):
        def __get__(self, instance, owner=None):
            return instance.remote_dict['value']

        def __set__(self, instance, val):
            instance.remote_dict['value'] = val

    value = IntDesc()

    def __init__(self, remote_dict, value=None):
        self.remote_dict = remote_dict
        if value is not None:
            self.value = value


class RemoteQueue(object):
    """Remote queue class, data saved in RemoteDict."""
    def __init__(self, remote_dict, value=None, mode='fifo'):
        self.remote_dict = remote_dict
        self.mode = mode
        if value is not None:
            self.remote_dict['value'] = value

    def get(self):
        val = self.remote_dict['value']
        if len(val) == 0:
            return None
        else:
            if self.mode == 'fifo':
                out = val[0]
                self.remote_dict['value'] = val[1:]
            elif self.mode == 'lifo':
                out = val[-1]
                self.remote_dict['value'] = val[:-1]

            return out

    def put(self, value):
        val = self.remote_dict['value']
        val.append(value)
        self.remote_dict['value'] = val


class RemoteDictMaster(dict):
    """Remote dictionary class - master side."""
    def set_name(self, name):
        self.name = name

    def remote_set(self, data, slave, mutable=False):
        key, value = data
        if not mutable and key in self:
            logger.error("imutable dict '%s'! key '%s' already in global dict"
                         % (self.name, key))
            mpi_comm.send(False, dest=slave, tag=tags.SET_DICT_STATUS)
            raise(KeyError)
        else:
            self.__setitem__(key, value)
            logger.debug('set master dict (%s[%s])' % (self.name, key))
            mpi_comm.send(True, dest=slave, tag=tags.SET_DICT_STATUS)

    def remote_get(self, key, slave):
        if key in self:
            mpi_comm.send(self.__getitem__(key), dest=slave, tag=tags.DICT_VAL)
            logger.debug('sent %s to %d (%s[%s])'
                         % (tags.name[tags.DICT_VAL], slave, self.name, key))
        else:
            mpi_comm.send(None, dest=slave, tag=tags.DICT_VAL)
            logger.error('RemoteDict KeyError (%s[%s])' % (self.name, key))
            raise(KeyError)

    def remote_get_keys(self, slave):
        mpi_comm.send(self.keys(), dest=slave, tag=tags.DICT_KEYS)
        logger.debug('sent %s to %d (%s)'
                     % (tags.name[tags.DICT_KEYS], slave, self.name))

    def remote_get_len(self, slave):
        mpi_comm.send(self.__len__(), dest=slave, tag=tags.DICT_LEN)
        logger.debug('sent %s to %d (%s)'
                     % (tags.name[tags.DICT_LEN], slave, self.name))

    def remote_get_in(self, key, slave):
        mpi_comm.send(self.__contains__(key), dest=slave, tag=tags.DICT_IN)
        logger.debug('sent %s to %d (%s)'
                     % (tags.name[tags.DICT_IN], slave, self.name))


class RemoteDict(object):
    """Remote dictionary class - slave side."""
    def __init__(self, name, mutable=False):
        self._dict = {}
        self.mutable = mutable
        self.name = name

    def _setitem(self, key, value, tag):
        mpi_comm.send((self.name, key, value), dest=mpi_master, tag=tag)
        logger.debug('sent %s to %d (%s[%s])'
                     % (tags.name[tag], mpi_master, self.name, key))
        stat = mpi_comm.recv(source=mpi_master, tag=tags.SET_DICT_STATUS)
        logger.debug('recevied %s from %d'
                     % (tags.name[tags.SET_DICT_STATUS], mpi_master))

        return stat

    def __setitem__(self, key, value):
        if self.mutable:
            self._setitem(key, value, tags.SET_DICT)
            self._dict[key] = value
        else:
            if key in self._dict:
                msg = "imutable dict '%s'! key '%s' already in local dict"\
                    % (self.name, key)
                logger.error(msg)
                raise(KeyError)
            else:
                stat = self._setitem(key, value, tags.SET_DICT_IMMUTABLE)
                if stat:
                    self._dict[key] = value
                else:
                    msg = \
                        "imutable dict '%s'! key '%s' already in global dict"\
                        % (self.name, key)
                    logger.error(msg)
                    raise(KeyError)

    def __getitem__(self, key):
        if key in self._dict and not self.mutable:
            logger.debug('get value from local dict! (%s[%s])'
                         % (self.name, key))
        else:
            mpi_comm.send((self.name, key), dest=mpi_master,
                          tag=tags.GET_DICT)
            logger.debug('sent %s to %d (%s)'
                         % (tags.name[tags.GET_DICT], mpi_master, key))
            data = mpi_comm.recv(source=mpi_master, tag=tags.DICT_VAL)
            logger.debug('received %s from %d' % (tags.name[tags.DICT_VAL],
                                                  mpi_master))
            if data is not None:
                self._dict[key] = data
            else:
                logger.error('RemoteDict KeyError (%s[%s])' % (self.name, key))
                raise(KeyError)

        return self._dict[key]

    def __len__(self):
        mpi_comm.send((self.name,), dest=mpi_master, tag=tags.GET_DICT_LEN)
        length = mpi_comm.recv(source=mpi_master, tag=tags.DICT_LEN)
        return length

    def __contains__(self, key):
        mpi_comm.send((self.name, key), dest=mpi_master, tag=tags.GET_DICT_IN)
        is_in = mpi_comm.recv(source=mpi_master, tag=tags.DICT_IN)
        return is_in

    def keys(self):
        mpi_comm.send((self.name,), dest=mpi_master, tag=tags.GET_DICT_KEYS)
        keys = mpi_comm.recv(source=mpi_master, tag=tags.DICT_KEYS)
        return keys

    def get(self, key, default=None):
        if key in self.keys():
            return self.__getitem__(key)
        else:
            return default

    def update(self, other):
        for k in other.keys():
            self.__setitem__(k, other[k])


class RemoteLock(object):
    """Remote lock class - lock and unlock restricted access to the master."""
    def __init__(self):
        self.locked = False

    def acquire(self):
        mpi_comm.send(None, dest=mpi_master, tag=tags.LOCK)
        data = mpi_comm.recv(source=mpi_master, tag=tags.LOCK_STATUS)
        if not data:
            logger.error('lock can not be aquired!')
            raise(SystemError)
        else:
            self.locked = True

    def release(self):
        mpi_comm.send(None, dest=mpi_master, tag=tags.UNLOCK)
        self.locked = False


def wait_for_tag(wtag, num=1):
    ndone = num
    while ndone > 0:
        mpi_comm.recv(source=MPI.ANY_SOURCE, tag=wtag, status=mpi_status)
        tag = mpi_status.Get_tag()
        source = mpi_status.Get_source()
        logger.debug('received %s from %d' % (tags.name[tag], source))
        if tag == wtag:
            ndone -= 1


def get_slaves():
    """Get the list of slave nodes"""
    slaves = range(mpi_comm.Get_size())
    slaves.remove(mpi_master)
    return slaves


def master_send_task(task, data, wait=False):
    """Send task to all slaves."""
    slaves = get_slaves()
    wait_for_tag(tags.READY, len(slaves))
    logger.info('all nodes are ready for task "%s"' % task)

    for ii in slaves:
        mpi_comm.send((task, data), dest=ii, tag=tags.START)

    if wait:
        wait_for_tag(tags.DONE, len(slaves))


def master_send_continue():
    """Send 'continue' to all slaves."""
    for ii in get_slaves():
        mpi_comm.send(None, dest=ii, tag=tags.CONTINUE)
    logger.info('slave nodes: continue')


def master_wait_for_done():
    """Wait for slaves to be done."""
    wait_for_tag(tags.DONE, len(get_slaves()))


def master_loop():
    """Run the master loop - wait for requests from slaves."""
    slaves = get_slaves()
    wait_for_tag(tags.READY, len(slaves))
    logger.info('all nodes are ready')
    for ii in slaves:
        mpi_comm.send(None, dest=ii, tag=tags.START)

    ndone = len(slaves)
    source = MPI.ANY_SOURCE
    while ndone > 0:
        data = mpi_comm.recv(source=source, tag=MPI.ANY_TAG, status=mpi_status)
        tag = mpi_status.Get_tag()
        slave = mpi_status.Get_source()
        logger.debug('received %s from %d' % (tags.name[tag], slave))
        if tag == tags.DONE:
            ndone -= 1
        elif tag == tags.LOCK:
            source = slave
            mpi_comm.send(True, dest=slave, tag=tags.LOCK_STATUS)
        elif tag == tags.UNLOCK:
            source = MPI.ANY_SOURCE
        elif tag == tags.SET_DICT:
            global_multiproc_dict[data[0]].remote_set(data[1:], slave,
                                                      mutable=True)
        elif tag == tags.SET_DICT_IMMUTABLE:
            global_multiproc_dict[data[0]].remote_set(data[1:], slave)
        elif tag == tags.GET_DICT:
            global_multiproc_dict[data[0]].remote_get(data[1], slave)
        elif tag == tags.GET_DICT_KEYS:
            global_multiproc_dict[data[0]].remote_get_keys(slave)
        elif tag == tags.GET_DICT_LEN:
            global_multiproc_dict[data[0]].remote_get_len(slave)
        elif tag == tags.GET_DICT_IN:
            global_multiproc_dict[data[0]].remote_get_in(data[1], slave)

    for ii in slaves:
        mpi_comm.send(None, dest=ii, tag=tags.CONTINUE)
    logger.info('slave nodes: continue')


def start_slave(tag=''):
    """Start the slave nodes."""
    mpi_comm.send(mpi_rank, dest=mpi_master, tag=tags.READY)
    logger.debug('%s ready' % tag)
    data = mpi_comm.recv(source=mpi_master, tag=tags.START)
    logger.info('%s started' % tag)

    return data


def stop_slave(tag='', wait=False):
    """Stop the slave nodes."""
    mpi_comm.send(mpi_rank, dest=mpi_master, tag=tags.DONE)
    logger.info('%s stoped' % tag)

    if wait:
        mpi_comm.recv(source=mpi_master, tag=tags.CONTINUE)
