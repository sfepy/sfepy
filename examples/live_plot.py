import sys
sys.path.append( '.' )

import time
from sfepy.base.base import *
from sfepy.base.log import Log


def main():

    log_conf = {
        'is_plot' : True,
    }

    log = Log.from_conf( log_conf, (['sin(x)'],) )

    for x in nm.linspace( 0, 4.0 * nm.pi, 200 ):
        log( nm.sin( x ) )
#        time.sleep( 0.05 )

    log( nm.sin( x ), finished = True )

    print log
    
if __name__ == '__main__':
    main()
