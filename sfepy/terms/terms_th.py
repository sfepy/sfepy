import numpy as nm

from sfepy.base.base import Struct
from sfepy.terms.terms import Term

class THTerm(Term):
    """
    Base class for terms depending on time history (fading memory
    terms).
    """

    def eval_real(self, shape, fargs, mode='eval', term_mode=None,
                  diff_var=None, **kwargs):

        if diff_var is None:
            if mode == 'eval':
                out = 0.0

            else:
                out = nm.zeros(shape, dtype=nm.float64)

            iter_kernel = fargs
            for ii, fargs in iter_kernel():
                out1, status = Term.eval_real(self, shape, fargs, mode=mode,
                                              term_mode=term_mode,
                                              diff_var=diff_var, **kwargs)
                out += out1

        else:
            out, status = Term.eval_real(self, shape, fargs, mode=mode,
                                         term_mode=term_mode,
                                         diff_var=diff_var, **kwargs)

        return out, status

class ETHTerm(Term):
    """
    Base class for terms depending on time history with exponential
    convolution kernel (fading memory terms).
    """

    def get_eth_data(self, key, state, decay, values):
        step_cache = state.evaluate_cache.setdefault('eth', {})
        cache = step_cache.setdefault(None, {})

        data_key = key + (self.arg_derivatives[state.name],)
        if data_key in cache:
            out = cache[data_key]
            out.values = values

        else:
            out = Struct(history=nm.zeros_like(values),
                         values=values,
                         decay=decay,
                         __advance__=self.advance_eth_data)
            cache[data_key] = out

        return out

    def advance_eth_data(self, ts, data):
        data.history[:] = data.decay * (data.history + data.values)
