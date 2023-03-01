from __future__ import absolute_import
import numpy as nm

from sfepy.base.base import assert_, OneTypeList, Container, Struct
import six

class Functions(Container):
    """Container to hold all user-defined functions."""

    def from_conf(conf):
        objs = OneTypeList(Function)
        for key, fc in six.iteritems(conf):
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

def make_sfepy_function(fun_or_name=None):
    """
    Convenience decorator to quickly create
    :class:`sfepy.discrete.functions.Function` objects.

    Has two modes of use either without parameter::

        @make_sfepy_function
        def my_function(...):
            ...

    or with name::

        @make_sfepy_function("new_name_for_my_function")
        def my_function(...):
            ...

    Parameters
    ----------
    fun_or_name : string, optional
        Name to be saved within `Function` instance,
        if None name of decorated function is used.

    Returns
    -------
    new_fun : `sfepy.discrete.functions.Function`
        With attribute name set to provided name or original function name.
    """
    if callable(fun_or_name):
        return Function(fun_or_name.__name__, fun_or_name)

    def functionizer(fun):
        """
        Internal decorator.

        Parameters
        ----------
        fun : callable
            Will be converted to sfepy.discrete.functions.Function.

        Returns
        -------
        new_fun : Function instance
            The `sfepy.discrete.functions.Function` object.
        """
        if fun_or_name is not None:
            return Function(fun_or_name, fun)
        return Function(fun.__name__, fun)
    return functionizer

class ConstantFunction(Function):
    """Function with constant values."""

    def __init__(self, values, no_tile=False):
        """Make a function out of a dictionary of constant values. When
        called with coors argument, the values are repeated for each
        coordinate. Not repeated if `no_tile` is True."""

        name = '_'.join(['get_constants'] + list(values.keys()))

        def get_constants(ts=None, coors=None, mode=None, **kwargs):
            out = {}
            if mode == 'special':
                for key, val in six.iteritems(values):
                    if '.' in key:
                        vkey = key.split('.')[1]
                        out[vkey] = val

            elif (mode == 'qp'):
                for key, val in six.iteritems(values):
                    if '.' in key: continue

                    dtype = nm.float64 if nm.isrealobj(val) else nm.complex128
                    val = nm.array(val, dtype=dtype, ndmin=3)
                    nc = 1 if no_tile else coors.shape[0]
                    out[key] = nm.tile(val, (nc, 1, 1))

            elif (mode == 'special_constant') or (mode is None):
                for key, val in six.iteritems(values):
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

        name = '_'.join(['get_constants_by_region'] + list(values.keys()))

        def get_constants(ts=None, coors=None, mode=None,
                          term=None, problem=None, **kwargs):
            out = {}
            if mode == 'qp':
                qps = term.get_physical_qps()
                assert_(qps.num == coors.shape[0])

                for key, val in six.iteritems(values):
                    if '.' in key: continue
                    rval = nm.array(val[list(val.keys())[0]], ndmin=3)
                    s0 = rval.shape[1:]
                    dtype = nm.float64 if nm.isrealobj(rval) else nm.complex128
                    matdata = nm.zeros(qps.shape[:2] + s0, dtype=dtype)

                    for rkey, rval in six.iteritems(val):
                        region = problem.domain.regions[rkey]
                        rval = nm.array(rval, dtype=dtype, ndmin=3)

                        cells = region.get_cells(true_cells_only=False)
                        tcells = term.region.get_cells(true_cells_only=False)
                        ii = nm.isin(tcells, cells)
                        matdata[ii] = rval

                    out[key] = matdata.reshape((-1,) + s0)

            return out

        Function.__init__(self, name=name, function=get_constants,
                          is_constant=True)
