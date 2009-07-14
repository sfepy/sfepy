from sfepy.base.base import *

class Functions(Container):
    """Container to hold all user-defined functions."""

    def from_conf(conf):
        objs = OneTypeList(Function)
        for key, fc in conf.iteritems():
            fun = Function(name = fc.name,
                           function = fc.function,
                           is_constant = False)
            objs.append(fun)

        obj = Functions(objs)
        return obj
    from_conf = staticmethod(from_conf)

class Function(Struct):
    """Base class for user-defined functions."""

    def __call__(self, *args, **kwargs):
        return self.function(*args, **kwargs)

class ConstantFunction(Function):

    def __init__(self, values, functions=None):
        """Make a function out of a dictionary of constant values."""

        name = '_'.join(['get_constants'] + values.keys())

        def get_constants(ts, coors, **kwargs):
            out = {}
            for key, val in values.iteritems():
                out[key] = val
            return out
        
        Function.__init__(self, name = name, function = get_constants,
                          is_constant = True)

        if functions is not None:
            functions.append(self)
