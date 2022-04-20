import os
import numpy as nm
import pytest

import sfepy.base.testing as tst

@pytest.fixture(scope='module')
def log_filename(output_dir):
    return os.path.join(output_dir, 'test_log.txt')

@pytest.fixture(scope='module')
def log(log_filename):
    from sfepy.base.log import Log

    log = Log([['x'], ['x^2', 'x^3']],
              plot_kwargs = [{},
                             [{'color' : 'b', 'ls' : '', 'marker' : 'o'},
                              {'color' : 'g', 'ls' : ':', 'marker' : 'x'}]],
              yscales=['log', 'linear'],
              xlabels=['x', 'x'], ylabels=['x', 'x^p'],
              is_plot=False,
              aggregate=0, sleep=0.0,
              log_filename=log_filename,
              formats=[['{:.3e}'], ['{:.5e}'] * 2])

    for x in nm.linspace(0, 1, 11):
        log(x, x**2, x**3, x=[x + 1, x])
        if nm.allclose(x, 0.5):
            log.plot_vlines([0], color='g', linewidth=2)
        if nm.allclose(x, 0.7):
            log.plot_vlines([1], color='g', linewidth=2)

    log(finished=True)
    return log

def test_log_rw(log_filename, log, output_dir):
    from sfepy.base.base import Output
    from sfepy.base.log import read_log, write_log

    log, info = read_log(log_filename)

    filename2 = os.path.join(output_dir, 'test_log2.txt')
    output = Output('', filename=filename2, quiet=True)

    write_log(output, log, info)

    log2, info2 = read_log(filename2)

    ok = True
    _ok = info == info2
    if not _ok:
        tst.report('headers are not equal!')
        tst.report(info)
        tst.report(info2)
    ok = ok and _ok

    for key, val2 in log2.items():
        val = log[key]
        _ok = nm.allclose(val[0], val2[0], rtol=0.0, atol=1e-14)
        if not _ok:
            tst.report('x values are not equal!')
            tst.report(val[0])
            tst.report(val2[0])
        ok = ok and _ok

        _ok = nm.allclose(val[1], val2[1], rtol=0.0, atol=1e-14)
        if not _ok:
            tst.report('y values are not equal!')
            tst.report(val[1])
            tst.report(val2[1])
        ok = ok and _ok

        _ok = nm.allclose(val[2], val2[2], rtol=0.0, atol=1e-14)
        if not _ok:
            tst.report('vlines are not equal!')
            tst.report(val[2])
            tst.report(val2[2])
        ok = ok and _ok

    assert ok
