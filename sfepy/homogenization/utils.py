from sfepy.base.base import *

##
# c: 28.02.2007, r: 13.02.2008
def build_op_pi( var_name, problem, ir, ic ):
    """\Pi^{rs}_i = \delta_{ri} y_s. """
    var = problem.variables[var_name]
    coor = var.field.get_coor()[:,:-1]

    pi = nm.zeros_like( coor )
    pi[:,ir] = coor[:,ic]
    pi.shape = (pi.shape[0] * pi.shape[1],)

    return pi

##
# c: 05.05.2008, r: 22.07.2008
def createPis( problem, variables, varName ):
    problem.setVariables( variables )

    dim = problem.domain.mesh.dim
    pis = nm.zeros( (dim, dim), dtype = nm.object )
    for ir in range( dim ):
        for ic in range( dim ):
            pi = buildOpPi( varName, problem, ir, ic )
            pis[ir,ic] = pi
    return pis

def iter_sym( dim ):
    for ii in xrange( dim ):
        yield ii, ii
    for ir in xrange( 0, dim ):
        for ic in xrange( ir + 1, dim ):
            yield ir, ic
c2s = {
    2 : [0, 2, 2, 1],
    3 : [0, 3, 4, 3, 1, 5, 4, 5, 2],
}
def coor_to_sym( ir, ic, dim ):
    return c2s[dim][dim*ir+ic]

##
# c: 14.09.2006, r: 04.04.2008
def interp_conv_mat( mat, ts, tdiff ):
    n_t = mat.shape[0]
    out = []
    tn = ts.time
    for ii, step in enumerate( xrange( ts.step, 0, -1 ) ):
        if ii == 0:
            out.append( mat[0] )
            continue
        
        td = tn - ts.times[step]
        if (td - 1e-12) > tdiff[-1]: break
        
        i1 = (tdiff >= td).argmax()
        i0 = i1 - 1

        td0, td1 = tdiff[[i0, i1]]
        dt = (td1 - td0)
        c1, c0 = (td - td0) / dt, (td1 - td) / dt
        out.append( c0 * mat[i0] + c1 * mat[i1] )

##         print ii, step, td
##        print i0, i1
##         print tn, ts.times[step]
##         print td0, td1, c0, c1

##     print out
##     print tdiff[-1], len( out )
##     pause()

    if not out: # For step == 0 matrix evaluation.
        out.append( mat[0] )

    return out

##
# c: 09.06.2008, r: 16.06.2008
def integrate_in_time( coef, ts, scheme = 'forward' ):
    """Forward difference or trapezoidal rule. 'ts' can be anything with
    'times' attribute."""
    if scheme == 'trapezoid':
        icoef = nm.sum( 0.5 * (coef[1:,...] + coef[:-1,...])
                        * nm.diff( ts.times )[:,nm.newaxis], axis = 0 )
    elif scheme == 'forward':
        icoef = nm.sum( coef[:-1,...]
                        * nm.diff( ts.times )[:,nm.newaxis], axis = 0 )
    else:
        raise ValueError( 'unsupported scheme: %s' % scheme )
    
    return icoef
