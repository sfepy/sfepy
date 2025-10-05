"""
Multiprocessing functions.
"""
import logging
import os

try:
    from mpi4py import MPI
    mpi_comm = MPI.COMM_WORLD
    mpi_rank = mpi_comm.Get_rank()
    mpi_status = MPI.Status()
    use_multiprocessing = mpi_comm.Get_size() > 1
except:
    use_multiprocessing = False


global_multiproc_dict = {}
mpi_master = 0


class MPILogFile:
    def __init__(self, comm, filename, mode):
        self._file = MPI.File.Open(comm, filename, mode)
        self._file.Set_atomicity(True)

    def write(self, msg):
        try:
            msg = msg.encode()
        except AttributeError:
            pass

        self._file.Write_shared(msg)

    def sync(self):
        self._file.Sync()

    def close(self):
        self.sync()
        self._file.Close()


class MPIFileHandler(logging.StreamHandler):
    "MPI file class for logging process communication."
    def __init__(self, filename, mode=MPI.MODE_WRONLY, comm=MPI.COMM_WORLD):
        self.baseFilename = os.path.abspath(filename)
        self.mode = mode
        self.comm = comm

        super(MPIFileHandler, self).__init__(self._open())

    def _open(self):
        stream = MPILogFile(self.comm, self.baseFilename, self.mode)
        return stream

    def close(self):
        if self.stream:
            self.stream.close()
            self.stream = None

    def emit(self, record):
        msg = self.format(record)
        self.stream.write('{}{}'.format(msg, self.terminator))
        self.flush()


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
                                  + ' - %(name)s: %(message)s')
    mh.setFormatter(formatter)
    logger.addHandler(mh)

    return logger


logger = get_logger()


def enum(*sequential):
    enums = dict(zip(sequential, range(len(sequential))))
    reverse = dict((value, key) for key, value in enums.items())
    enums['name'] = reverse
    return type('Enum', (), enums)


tags = enum('READY', 'DONE', 'START', 'CONTINUE',
            'LOCK', 'UNLOCK',
            'SET_DICT', 'SET_DICT_IMMUTABLE', 'GET_DICT', 'DICT_VAL',
            'SET_DICT_STATUS', 'GET_DICT_KEYS', 'DICT_KEYS',
            'GET_DICT_LEN', 'DICT_LEN', 'GET_DICT_IN', 'DICT_IN',
            'GET_QUEUE', 'PUT_QUEUE', 'QUEUE_VAL')


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
    if mpi_rank == mpi_master:
        if key not in global_multiproc_dict:
            queue = RemoteQueueMaster(name)
            global_multiproc_dict[key] = queue
        else:
            queue = global_multiproc_dict[key]
            queue.clean()
    else:
        queue = RemoteQueue(name)

    return queue


def get_dict(name, mutable=False, clear=False, soft_set=False):
    """Get the remote dictionary."""
    if mpi_rank == mpi_master:
        if name in global_multiproc_dict:
            if clear:
                global_multiproc_dict[name].clear()
                logger.info("cleaning dict %s" % name)
            else:
                logger.info("using existing dict %s" % name)
            return global_multiproc_dict[name]
        else:
            rdict = RemoteDictMaster(name, mutable=mutable, soft_set=soft_set)
            global_multiproc_dict[name] = rdict
            logger.info("new dict %s" % name)
            return rdict
    else:
        return RemoteDict(name, mutable=mutable)


class RemoteInt:
    """Remote intiger class, data saved in RemoteDict."""
    class IntDesc:
        def __get__(self, instance, owner=None):
            return instance.remote_dict['value']

        def __set__(self, instance, val):
            instance.remote_dict['value'] = val

    value = IntDesc()

    def __init__(self, remote_dict, value=None):
        self.remote_dict = remote_dict
        if value is not None:
            self.value = value


class RemoteQueueMaster(list):
    """Remote queue class - master side."""
    def __init__(self, name, mode='fifo', *args):
        list.__init__(self, *args)
        self.name = name
        self.mode = mode

    @staticmethod
    def get_gdict_key(name):
        return '__queue_%s' % name

    def get(self):
        if len(self) == 0:
            return None
        else:
            if self.mode == 'fifo':
                out = self.pop(0)
            elif self.mode == 'lifo':
                out = self.pop(-1)

            return out

    def put(self, value):
        self.append(value)

    def clean(self):
        del self[:]

    def remote_put(self, value, slave):
        self.put(value)

    def remote_get(self, slave):
        mpi_comm.isend(self.get(), dest=slave, tag=tags.QUEUE_VAL)
        logger.debug('sent %s to %d (%s)'
                     % (tags.name[tags.QUEUE_VAL], slave, self.name))


class RemoteQueue:
    """Remote queue class - slave side."""
    def __init__(self, name):
        self.name = name

    def get(self):
        mpi_comm.isend(self.name, dest=mpi_master, tag=tags.GET_QUEUE)
        logger.debug('sent %s to %d (%s)'
                     % (tags.name[tags.GET_QUEUE], mpi_master, self.name))
        value = mpi_comm.recv(source=mpi_master, tag=tags.QUEUE_VAL)
        logger.debug('received %s from %d'
                     % (tags.name[tags.QUEUE_VAL], mpi_master))

        return value

    def put(self, value):
        mpi_comm.isend((self.name, value), dest=mpi_master, tag=tags.PUT_QUEUE)
        logger.debug('sent %s to %d (%s)'
                     % (tags.name[tags.PUT_QUEUE], mpi_master, self.name))


class RemoteDictMaster(dict):
    """Remote dictionary class - master side."""
    def __init__(self, name, mutable=False, soft_set=False, *args):
        dict.__init__(self, *args)
        self.name = name
        self.mutable = mutable
        self.immutable_soft_set = soft_set

    def remote_set(self, data, slave, mutable=False):
        key, value = data
        if not self.mutable and key in self:
            if self.immutable_soft_set:
                mpi_comm.isend(False, dest=slave, tag=tags.SET_DICT_STATUS)
            else:
                msg = "imutable dict '%s'! key '%s' already in global dict"\
                    % (self.name, key)
                logger.error(msg)
                mpi_comm.isend(False, dest=slave, tag=tags.SET_DICT_STATUS)
                raise(KeyError)
        else:
            self.__setitem__(key, value)
            logger.debug('set master dict (%s[%s])' % (self.name, key))
            mpi_comm.isend(True, dest=slave, tag=tags.SET_DICT_STATUS)

    def remote_get(self, key, slave):
        if key in self:
            mpi_comm.isend(self.__getitem__(key), dest=slave,
                           tag=tags.DICT_VAL)
            logger.debug('sent %s to %d (%s[%s])'
                         % (tags.name[tags.DICT_VAL], slave, self.name, key))
        else:
            mpi_comm.isend(None, dest=slave, tag=tags.DICT_VAL)
            logger.error('RemoteDict KeyError (%s[%s])' % (self.name, key))
            raise(KeyError)

    def remote_get_keys(self, slave):
        mpi_comm.isend(self.keys(), dest=slave, tag=tags.DICT_KEYS)
        logger.debug('sent %s to %d (%s)'
                     % (tags.name[tags.DICT_KEYS], slave, self.name))

    def remote_get_len(self, slave):
        mpi_comm.isend(self.__len__(), dest=slave, tag=tags.DICT_LEN)
        logger.debug('sent %s to %d (%s)'
                     % (tags.name[tags.DICT_LEN], slave, self.name))

    def remote_get_in(self, key, slave):
        mpi_comm.isend(self.__contains__(key), dest=slave, tag=tags.DICT_IN)
        logger.debug('sent %s to %d (%s)'
                     % (tags.name[tags.DICT_IN], slave, self.name))


class RemoteDict:
    """Remote dictionary class - slave side."""
    def __init__(self, name, mutable=False):
        self._dict = {}
        self.mutable = mutable
        self.name = name

    def _setitem(self, key, value, tag):
        mpi_comm.isend((self.name, key, value), dest=mpi_master, tag=tag)
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
            mpi_comm.isend((self.name, key), dest=mpi_master,
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
        mpi_comm.isend((self.name,), dest=mpi_master, tag=tags.GET_DICT_LEN)
        length = mpi_comm.recv(source=mpi_master, tag=tags.DICT_LEN)
        return length

    def __contains__(self, key):
        if key in self._dict and not self.mutable:
            is_in = True
        else:
            mpi_comm.isend((self.name, key), dest=mpi_master,
                           tag=tags.GET_DICT_IN)
            is_in = mpi_comm.recv(source=mpi_master, tag=tags.DICT_IN)

        return is_in

    def keys(self):
        mpi_comm.isend((self.name,), dest=mpi_master, tag=tags.GET_DICT_KEYS)
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


def is_remote_dict(d):
    """Return True if 'd' is RemoteDict or RemoteDictMaster instance."""
    return isinstance(d, RemoteDict) or isinstance(d, RemoteDictMaster)


class RemoteLock:
    """Remote lock class - lock and unlock restricted access to the master."""
    def __init__(self):
        self.locked = False

    def acquire(self):
        mpi_comm.isend(None, dest=mpi_master, tag=tags.LOCK)
        self.locked = True

    def release(self):
        mpi_comm.isend(None, dest=mpi_master, tag=tags.UNLOCK)
        self.locked = False


def wait_for_tag(wtag, num=1):
    ndone = num
    start = MPI.Wtime()
    while ndone > 0:
        mpi_comm.recv(source=MPI.ANY_SOURCE, tag=wtag, status=mpi_status)
        tag = mpi_status.Get_tag()
        source = mpi_status.Get_source()
        logger.debug('received %s from %d (%.03fs)' % (tags.name[tag],
                                                       source,
                                                       MPI.Wtime() - start))
        if tag == wtag:
            ndone -= 1


def slave_get_task(name=''):
    """Start the slave nodes."""
    mpi_comm.isend(mpi_rank, dest=mpi_master, tag=tags.READY)
    logger.debug('%s ready' % name)
    task, data = mpi_comm.bcast(None, root=mpi_master)
    logger.info('%s received task %s' % (name, task))

    return task, data


def slave_task_done(task=''):
    """Stop the slave nodes."""
    mpi_comm.isend(mpi_rank, dest=mpi_master, tag=tags.DONE)
    logger.info('%s stopped' % task)


def get_slaves():
    """Get the list of slave nodes"""
    slaves = list(range(mpi_comm.Get_size()))
    slaves.remove(mpi_master)
    return slaves


def master_send_task(task, data):
    """Send task to all slaves."""
    slaves = get_slaves()
    wait_for_tag(tags.READY, len(slaves))
    logger.info('all nodes are ready for task %s' % task)
    mpi_comm.bcast((task, data), root=mpi_master)


def master_send_continue():
    """Send 'continue' to all slaves."""
    for ii in get_slaves():
        mpi_comm.send(None, dest=ii, tag=tags.CONTINUE)
    logger.info('slave nodes - continue')


def master_loop():
    """Run the master loop - wait for requests from slaves."""
    logger.info('main loop started')
    master_send_task('calculate', None)

    ndone = len(get_slaves())
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
        elif tag == tags.GET_QUEUE:
            qkey = RemoteQueueMaster.get_gdict_key(data)
            global_multiproc_dict[qkey].remote_get(slave)
        elif tag == tags.PUT_QUEUE:
            qkey = RemoteQueueMaster.get_gdict_key(data[0])
            global_multiproc_dict[qkey].remote_put(data[1], slave)

    logger.info('main loop finished')
