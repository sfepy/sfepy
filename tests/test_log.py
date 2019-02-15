import os
import numpy as nm

from sfepy.base.testing import TestCommon

class Test(TestCommon):
    tests = ['test_log_create', 'test_log_rw']

    @staticmethod
    def from_conf(conf, options):
        return Test(
            conf=conf, options=options,
            log_filename=os.path.join(options.out_dir, 'test_log.txt'),
        )

    def test_log_create(self):
        from sfepy.base.log import Log

        log = Log([['x^3']],
                  plot_kwargs = [[{'color' : 'b', 'ls' : '', 'marker' : 'o'}]],
                  yscales=['linear'],
                  xlabels=['x'], ylabels=['a cubic function'],
                  is_plot=False,
                  aggregate=0, sleep=0.0,
                  log_filename=self.log_filename,
                  formats=[['{:.5e}']])

        for x in nm.linspace(0, 1, 11):
            log(x**3, x=[x])
            if nm.allclose(x, 0.5):
                log.plot_vlines([0], color='g', linewidth=2)

        log(finished=True)

        return True

    def test_log_rw(self):
        from sfepy.base.base import Output
        from sfepy.base.log import read_log, write_log

        log, info = read_log(self.log_filename)

        output = Output('', filename=os.path.join(self.options.out_dir,
                                                  'test_log2.txt'),
                        quiet=True)
        write_log(output, log, info)

        log2, info2 = read_log(self.log_filename)

        ok = True
        _ok = info == info2
        if not _ok:
            self.report('headers are not equal!')
            self.report(info)
            self.report(info2)
        ok = ok and _ok

        for key, val2 in log2.items():
            val = log[key]
            _ok = nm.allclose(val[0], val2[0], rtol=0.0, atol=1e-14)
            if not _ok:
                self.report('x values are not equal!')
                self.report(val[0])
                self.report(val2[0])
            ok = ok and _ok

            _ok = nm.allclose(val[1], val2[1], rtol=0.0, atol=1e-14)
            if not _ok:
                self.report('y values are not equal!')
                self.report(val[1])
                self.report(val2[1])
            ok = ok and _ok

            _ok = nm.allclose(val[2], val2[2], rtol=0.0, atol=1e-14)
            if not _ok:
                self.report('vlines are not equal!')
                self.report(val[2])
                self.report(val2[2])
            ok = ok and _ok

        return ok
