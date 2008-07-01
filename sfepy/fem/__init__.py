try:
    from extmods import *
except:
    from sfepy.base.base import output
    output( 'warning: sfepy extension modules are not compiled!' )
    output( 'type "make"' )
