import numpy as nm

from sfepy.base.base import OneTypeList, Container, Struct

class Functions(Container):
    """Container to hold all user-defined functions."""

    def from_conf(conf):
        objs = OneTypeList(Function)
        for key, fc in conf.iteritems():
            fun = Function(name = fc.name,
                           function = fc.function,
                           is_constant = False,
                           extra_args = {})
            objs.append(fun)

        obj = Functions(objs)
        return obj
    from_conf = staticmethod(from_conf)

class Function(Struct):
    """Base class for user-defined functions."""

    def __init__(self, name, function, is_constant=False, extra_args=None):
        Struct.__init__(self, name = name, function = function,
                        is_constant = is_constant)
        if extra_args is None:
            extra_args = {}
        self.extra_args = extra_args

    def __call__(self, *args, **kwargs):
        _kwargs = dict(kwargs)
        _kwargs.update(self.extra_args)
        return self.function(*args, **_kwargs)

    def set_function(self, function, is_constant=False):
        self.function = function
        self.is_constant = is_constant

    def set_extra_args(self, **extra_args):
        self.extra_args = extra_args

class ConstantFunction(Function):
    """Function with constant values."""

    def __init__(self, values):
        """Make a function out of a dictionary of constant values. When
        called with coors argument, the values are repeated for each
        coordinate."""

        name = '_'.join(['get_constants'] + values.keys())

        def get_constants(ts=None, coors=None, mode=None, **kwargs):
            out = {}
            if mode == 'special':
                for key, val in values.iteritems():
                    if '.' in key:
                        vkey = key.split('.')[1]
                        out[vkey] = val

            elif (mode == 'qp'):
                for key, val in values.iteritems():
                    if '.' in key: continue

                    val = nm.array(val, dtype=nm.float64, ndmin=3)
                    out[key] = nm.tile(val, (coors.shape[0], 1, 1))

            elif (mode == 'special_constant') or (mode is None):
                for key, val in values.iteritems():
                    if '.' in key: continue

                    out[key] = val

            else:
                raise ValueError('unknown function mode! (%s)' % mode)
            return out

        Function.__init__(self, name = name, function = get_constants,
                          is_constant = True)

class ConstantFunctionByRegion(Function):
    """
    Function with constant values in regions.
    """

    def __init__(self, values):
        """
        Make a function out of a dictionary of constant values per region. When
        called with coors argument, the values are repeated for each
        coordinate in each of the given regions.
        """

        name = '_'.join(['get_constants_by_region'] + values.keys())

        def get_constants(ts=None, coors=None, mode=None,
                          term=None, problem=None, **kwargs):
            out = {}
            if mode == 'qp':
                qps = term.get_physical_qps()

                for key, val in values.iteritems():
                    if '.' in key: continue
                    rval = nm.array(val[val.keys()[0]], dtype=nm.float64,
                                    ndmin=3)
                    matdata = nm.zeros((coors.shape[0], ) + rval.shape[1:],
                                       dtype=nm.float64)

                    for rkey, rval in val.iteritems():
                        region = problem.domain.regions[rkey]
                        rval = nm.array(rval, dtype=nm.float64, ndmin=3)

                        for kgrp, elems in region.cells.iteritems():
                            nqp = qps.shape[kgrp][1]
                            nel = elems.shape[0]
                            vmap = nm.tile(elems.reshape((nel,1)) * nqp,
                                           (1, nqp)) + nm.arange(nqp)
                            matdata[vmap.reshape(nel * nqp)
                                    + qps.rindx[kgrp].start] = rval

                    out[key] = matdata

            return out

        Function.__init__(self, name=name, function=get_constants,
                          is_constant=True)
