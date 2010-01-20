import sys
sys.path.append( '.' )

from sfepy.base.base import *
from sfepy.base.log import Log

def main():
    log = Log((['sin( x )', 'cos( x )'], ['exp( x )']),
              yscales=['linear', 'log'],
              xlabels=['angle', None], ylabels=[None, 'a function'],
              log_filename='examples/live_plot.log')

    added = 0
    for x in nm.linspace( 0, 4.0 * nm.pi, 200 ):
        output( 'x: ', x )

        if x < (2.0 * nm.pi):
            log( nm.sin( x ), nm.cos( x ), nm.exp( x ), x = [x, None] )

        else:
            if added:
                log(nm.sin( x ), nm.cos( x ), nm.exp( x ), x**2,
                    x=[x, None, x])
            else:
                log.plot_vlines(color='r', linewidth=2)
                log.add_group(['x^2'], 'linear', 'new x', 'square',
                              formats=['%+g'])
            added += 1

        if (added == 20) or (added == 50):
            log.plot_vlines([2], color='g', linewidth=2)

    print log
    pause()

    log( finished = True )

if __name__ == '__main__':
    main()
