import sys
sys.path.append( '.' )

from sfepy.base.base import *
from sfepy.base.log import Log

def main():
    log = Log((['sin( x )', 'cos( x )'], ['exp( x )']),
              yscales=['linear', 'log'],
              xaxes=['angle', None], yaxes=[None, 'a function'])

    for x in nm.linspace( 0, 4.0 * nm.pi, 200 ):
        output( 'x: ', x )
        log( nm.sin( x ), nm.cos( x ), nm.exp( x ), x = [x, None] )

    print log
    pause()

    log( finished = True )

if __name__ == '__main__':
    main()
