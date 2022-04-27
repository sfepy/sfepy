from __future__ import print_function
from __future__ import absolute_import
import os
import sys
sys.path.append( '.' )

import numpy as nm

from sfepy.base.base import output, pause
from sfepy.base.log import Log

def main():
    cwd = os.path.split(os.path.join(os.getcwd(), __file__))[0]

    log = Log((['sin(x) + i sin(x**2)', 'cos(x)'], ['exp(x)']),
              yscales=['linear', 'log'],
              xlabels=['angle', None], ylabels=[None, 'a function'],
              aggregate=1000, sleep=0.05,
              log_filename=os.path.join(cwd, 'live_plot.log'))
    log2 = Log([['x^3']],
               yscales=['linear'],
               xlabels=['x'], ylabels=['a cubic function'],
               aggregate=1000, sleep=0.05,
               log_filename=os.path.join(cwd, 'live_plot2.log'),
               formats=[['{:.5e}']])

    added = 0
    for x in nm.linspace(0, 4.0 * nm.pi, 200):
        output('x: ', x)

        if x < (2.0 * nm.pi):
            log(nm.sin(x)+1j*nm.sin(x**2), nm.cos(x), nm.exp(x), x=[x, None])

        else:
            if added:
                log(nm.sin(x)+1j*nm.sin(x**2), nm.cos(x), nm.exp(x), x**2,
                    x=[x, None, x])
            else:
                log.plot_vlines(color='r', linewidth=2)
                log.add_group(['x^2'], yscale='linear', xlabel='new x',
                              ylabel='square', formats=['%+g'])
            added += 1

        if (added == 20) or (added == 50):
            log.plot_vlines([2], color='g', linewidth=2)

        log2(x*x*x, x=[x])

    print(log)
    print(log2)
    pause()

    log(finished=True)
    log2(finished=True)

if __name__ == '__main__':
    main()
