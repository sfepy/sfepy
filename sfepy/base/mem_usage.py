"""
Memory usage functions.
"""
from __future__ import absolute_import
import sys
import collections

import numpy as nm
import scipy.sparse as sp

from sfepy.base.base import basestr, Struct, Output
import six

def get_mem_usage(obj, usage=None, name=None, traversal_order=None, level=0):
    """
    Get lower bound of memory usage of an object.

    Takes into account strings, numpy arrays and scipy CSR sparse matrices,
    descends into sequences, mappings and objects.

    Parameters
    ----------
    obj : any object
        The object to be measured.
    usage : dict
        The dict with memory usage records, serving also as a cache of already
        traversed objects.
    name : str
        The name to be given to the object in its record.
    traversal_order : list, internal
        The traversal order of the object.
    level : int, internal
        The recurrence level.

    Returns
    -------
    usage : int
        The object's lower bound of memory usage.
    """
    if usage is None:
        usage = {}

    if name is None:
        name = getattr(obj, 'name', '-')

    if traversal_order is None:
        traversal_order = [0]

    to = traversal_order

    key = id(obj)
    if key in usage:
        usage[key].nrefs += 1
        return 0

    else:
        record = usage.setdefault(key, Struct(name=name,
                                              kind=type(obj).__name__,
                                              usage=0, nrefs=1,
                                              traversal_order=to[0],
                                              level=level))
        level += 1

    if isinstance(obj, nm.ndarray):
        record.usage = obj.nbytes

    elif isinstance(obj, sp.csr_matrix):
        record.usage = (get_mem_usage(obj.data, usage, name='data',
                                      traversal_order=to, level=level)
                        + get_mem_usage(obj.indices, usage, name='indices',
                                        traversal_order=to, level=level)
                        + get_mem_usage(obj.indptr, usage, name='indptr',
                                        traversal_order=to, level=level))

    elif isinstance(obj, basestr):
        record.usage = len(obj)

    elif isinstance(obj, Struct):
        for subname, sub in six.iteritems(obj.__dict__):
            to[0] += 1
            record.usage += get_mem_usage(sub, usage,
                                          name='attribute %s of %s'
                                          % (subname, getattr(obj, 'name',
                                                              record.kind)),
                                          traversal_order=to, level=level)

    elif isinstance(obj, collections.Mapping):
        try:
            for subname, sub in six.iteritems(obj):
                to[0] += 1
                record.usage += get_mem_usage(sub, usage,
                                              name='item %s of %s'
                                              % (subname, record.kind),
                                              traversal_order=to, level=level)
        except:
            pass

    elif isinstance(obj, collections.Sequence):
        for ii, sub in enumerate(obj):
            to[0] += 1
            record.usage += get_mem_usage(sub, usage,
                                          name='item %d of %s'
                                          % (ii, record.kind),
                                          traversal_order=to, level=level)

    else:
        record.usage = sys.getsizeof(obj)

    return record.usage

def print_mem_usage(usage, order_by='usage', direction='up', print_key=False):
    """
    Print memory usage dictionary.

    Parameters
    ----------
    usage : dict
        The dict with memory usage records.
    order_by : 'usage', 'name', 'kind', 'nrefs', 'traversal_order', or 'level'
        The sorting field name.
    direction : 'up' or 'down'
        The sorting direction.
    print_key : bool
        If True, print also the record key (object's id).
    """
    keys = list(usage.keys())
    order_vals = nm.array([record.get(order_by)
                           for record in six.itervalues(usage)])

    order = nm.argsort(order_vals)
    if direction == 'down':
        order = order[::-1]

    output = Output('')
    fmt = '%9s, %s, %s, %d %d %d' + ', %d' * print_key

    for ii in order:
        key = keys[ii]
        record = usage[key]

        if print_key:
            output(fmt % (record.usage, record.name, record.kind, record.nrefs,
                          record.traversal_order, record.level, key))

        else:
            output(fmt % (record.usage, record.name, record.kind, record.nrefs,
                          record.traversal_order, record.level))

def raise_if_too_large(size, factor=1.0):
    """
    Raise MemoryError if the total system memory is lower than `size` times
    safety `factor`. Use `factor=None` for skipping the memory check.
    """
    if factor is None: return

    import psutil
    mem = psutil.virtual_memory()
    if factor * size > mem.total:
        mb = 1000**2
        raise MemoryError('insufficent memory {} MB to allocate {} MB'
                          ' with safety factor {}'
                          .format(mem.total/mb, size/mb, factor))
