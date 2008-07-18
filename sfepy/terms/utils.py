from sfepy.base.base import *

##
# c: 31.07.2007, r: 27.03.2008
def fix_scalar_in_el( mat, n_el, dtype, default = 1.0 ):
##     print '>>>', mat
    if mat is None:
        out = nm.empty( (n_el, 1, 1, 1), dtype = dtype )
        out.fill( default )
    elif nm.isscalar( mat ):
        out = nm.empty( (n_el, 1, 1, 1), dtype = dtype )
        out.fill( mat )
    elif isinstance( mat, nm.ndarray ):
        if mat.size == n_el:
            out = mat.reshape( n_el, 1, 1, 1 )
        else:
            out = nm.empty( (n_el, 1, 1, 1), dtype = dtype )
            out.fill( mat )
    else:
        print 'unknown mat type!'
        print mat
        raise ValueError
        
    return out

##
# c: 27.03.2008, r: 27.03.2008
def fix_scalar_constant( mat, dtype ):
    out = None
    if nm.isscalar( mat ):
        out = dtype( mat )
    elif isinstance( mat, nm.ndarray ):
        if mat.size == 1:
            out = dtype( mat )
    return out

##
# c: 19.12.2007, r: 01.04.2008
def choose_scalar_or_in_el( mat, dtype, fun1, fun2 ):
    if nm.isscalar( mat ):
        out = dtype( mat ), fun1
    elif isinstance( mat, nm.ndarray ):
        if mat.ndim == 0:
            out = dtype( mat ), fun1
        else:
            out = mat, fun2
    return out

