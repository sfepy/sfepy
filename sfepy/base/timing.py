"""
Elapsed time measurement utilities.
"""
import time

from sfepy.base.base import Struct

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
