"""
Elapsed time measurement utilities.
"""
import time

from sfepy.base.base import PY3, Struct

class Timer(Struct):

    def __init__(self, name):
        Struct.__init__(self, name=name)
        self.time_function = time.perf_counter if PY3 else time.clock
        self.reset()

    def reset(self):
        self.t0 = self.t1 = None
        self.total = self.dt = 0.0

    def start(self, reset=False):
        if reset: self.reset()

        self.t0 = self.time_function()
        self.t1 = None

    def stop(self):
        if self.t0 is None:
            raise ValueError('timer %s was not started!' % self.name)
        self.t1 = self.time_function()

        self.dt = self.t1 - self.t0
        self.total += self.dt
        return self.dt
