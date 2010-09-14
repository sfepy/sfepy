from sfepy.terms.extmods import terms
from sfepy.terms.cache import DataCache
from sfepy.base.base import nm, pause, debug

class ExpHistoryDataCache(DataCache):
    """History for exponential decay convolution kernels.

    The decay argument is F(\Delta t), F(t_0=0) is assumed to be 1.0.
    """
    name = 'exp_history'
    arg_types = ('decay', 'values')

    def __init__(self, name, arg_names, history_sizes=None):
        DataCache.__init__(self, name, arg_names,
                           ['history', 'increment', 'decay'], history_sizes)

    def init_data(self, key, ckey, term, **kwargs):
        decay, values = self.get_args(**kwargs)
        shape = values.shape
        self.shapes = {
            'history' : shape,
            'increment' : shape,
            'decay' : decay.shape,
        }
        DataCache.init_datas(self, ckey, self.shapes, zero=True)

    def update(self, key, term, ih, **kwargs):
        decay, values = self.get_args(**kwargs)

        ckey = self.get_key(term)

        self.data['increment'][ckey][ih] = values
        self.data['decay'][ckey][ih] = decay

        self.valid['history'][ckey] = True
        self.valid['increment'][ckey] = True
        self.valid['decay'][ckey] = True

    def custom_advance(self, key, ckey, step):
        if key == 'history':
            history = self.data['history'][ckey][0]
            increment = self.data['increment'][ckey][0]
            decay = self.data['decay'][ckey][0]

            self.data['history'][ckey][0][:] = decay * (history + increment)
