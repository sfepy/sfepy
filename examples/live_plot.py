import sys
sys.path.append( '.' )

from sfepy.base.base import *
from sfepy.base.log import Log


def main():

    log_conf = {
        'is_plot' : True,
        'aggregate' : 200,
        'yscales' : ['linear', 'log'],
        'xaxes' : ['angle', None],
    }

    log = Log.from_conf( log_conf, (['sin( x )', 'cos( x )'],['exp( x )']) )

    for x in nm.linspace( 0, 4.0 * nm.pi, 200 ):
        output( 'x: ', x )
        log( nm.sin( x ), nm.cos( x ), nm.exp( x ), x = [x, None] )

    print log
    pause()

    log( finished = True )


if __name__ == '__main__':
    main()
