class TimeSuite:
    """
    An example of parametrized time benchmark.
    """
    params = [0, 10, 20, 1000]
    param_names = ['n']

    def setup(self, n):
        self.obj = range(n)

    def teardown(self, n):
        del self.obj

    def time_range_iter(self, n):
        for i in self.obj:
            print(i)
            pass


class MemSuite:
    """
    An example of memory benchmark.
    """
    def mem_list(self):
        return [0] * 256
