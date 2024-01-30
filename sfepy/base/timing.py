"""
Elapsed time measurement utilities.
"""
import time

from sfepy.base.base import Struct, IndexedStruct

class Timer(Struct):

    def __init__(self, name='timer', start=False):
        Struct.__init__(self, name=name)
        self.time_function = time.perf_counter
        self.reset()

        if start:
            self.start()

    def reset(self):
        self.t0 = self.t1 = None
        self.total = self.dt = 0.0

    def start(self, reset=False):
        if reset: self.reset()

        self.t1 = None
        self.t0 = self.time_function()

    def stop(self):
        self.t1 = self.time_function()
        if self.t0 is None:
            raise ValueError('timer "%s" was not started!' % self.name)

        self.dt = self.t1 - self.t0
        self.total += self.dt
        return self.dt

    def add(self, dt):
        if self.t1 is None:
            raise ValueError('timer "%s" was not started and stopped!'
                             % self.name)

        self.t1 += dt
        self.dt += dt
        self.total += dt

class Timers(IndexedStruct):

    def __init__(self, names):
        IndexedStruct.__init__(self,
                               **{name : Timer(name=name) for name in names})

    def create(self, name):
        if name in self.__dict__:
            self[name].reset()

        else:
            self[name] = Timer(name=name)

    def start(self, name):
        self[name].start()

    def stop(self, name):
        return self[name].stop()

    def get_dts(self):
        return {name : timer.dt for name, timer in self.to_dict().items()}

    def get_totals(self):
        return {name : timer.total for name, timer in self.to_dict().items()}
