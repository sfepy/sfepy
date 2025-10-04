import numpy as nm

def build_op_pi(var, ir, ic):
    r"""\Pi_i^{rs} = y_s \delta_{ir} for r = `ir`, s = `ic`."""
    coor = var.field.get_coor()

    pi = nm.zeros_like( coor )
    pi[:,ir] = coor[:,ic]
    pi.shape = (pi.shape[0] * pi.shape[1],)

    return pi

def create_pis(problem, var_name):
    r"""\Pi_i^{rs} = y_s \delta_{ir}, \ul{y} \in Y coordinates."""
    var = problem.get_variables(auto_create=True)[var_name]

    dim = problem.domain.mesh.dim
    pis = nm.zeros((dim, dim), dtype=object)
    components = []
    for ir in range( dim ):
        for ic in range( dim ):
            pi = build_op_pi(var, ir, ic)
            pis[ir,ic] = {var_name : pi}
            components.append((ir, ic))
    return components, pis

def create_scalar_pis( problem, var_name ):
    r"""\Pi^k = y_k, \ul{y} \in Y coordinates."""
    var = problem.get_variables(auto_create=True)[var_name]
    coor = var.field.get_coor()

    dim = problem.domain.mesh.dim
    pis = nm.zeros((dim,), dtype=object)
    components = []
    for ir in range( dim ):
        pis[ir] = {var_name : nm.ascontiguousarray( coor[:,ir] )}
        components.append((ir,))
    return components, pis

def iter_sym( dim ):
    for ii in range( dim ):
        yield ii, ii
    for ir in range( 0, dim ):
        for ic in range( ir + 1, dim ):
            yield ir, ic

def iter_nonsym(dim):
    for ir in range(dim):
        for ic in range(dim):
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
    for ii, step in enumerate( range( ts.step, 0, -1 ) ):
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
    dt = nm.diff(ts.times)
    dt = dt.reshape((dt.shape[0],) + (1,) * (coef.ndim-1))

    if scheme == 'trapezoid':
        icoef = nm.sum(0.5 * (coef[1:,...] + coef[:-1,...]) * dt, axis=0)
    elif scheme == 'forward':
        icoef = nm.sum(coef[:-1,...] * dt, axis=0)
    else:
        raise ValueError( 'unsupported scheme: %s' % scheme )

    return icoef

def define_box_regions(dim, lbn, rtf=None, eps=1.0e-3, kind='facet'):
    """
    Define sides and corner regions for a box aligned with coordinate
    axes.

    Parameters
    ----------
    dim : int
        Space dimension
    lbn : tuple
        Left bottom near point coordinates if rtf is not None. If rtf is
        None, lbn are the (positive) distances from the origin.
    rtf : tuple
        Right top far point coordinates.
    eps : float
        A parameter, that should be smaller than the smallest mesh node
        distance.
    kind : bool, optional
       The region kind.

    Returns
    -------
    regions : dict
        The box regions.
    """
    if rtf is None:
        lbn, rtf = -nm.array(lbn), lbn

    if dim == 3:
        lbnx, lbny, lbnz = lbn
        rtfx, rtfy, rtfz = rtf
        dx = abs(rtfx-lbnx)
        dy = abs(rtfy-lbny)
        dz = abs(rtfz-lbnz)
        lbnx, lbny, lbnz = (lbnx+dx*eps, lbny+dy*eps, lbnz+dz*eps)
        rtfx, rtfy, rtfz = (rtfx-dx*eps, rtfy-dy*eps, rtfz-dz*eps)
        regions = {
            'Near' : ('vertices in (y < %.16e)' % lbny, kind),
            'Far' : ('vertices in (y > %.16e)' % rtfy, kind),
            'Bottom' : ('vertices in (z < %.16e)' % lbnz, kind),
            'Top' : ('vertices in (z > %.16e)' % rtfz, kind),
            'Left' : ('vertices in (x < %.16e)' % lbnx, kind),
            'Right' : ('vertices in (x > %.16e)' % rtfx, kind),
            'Corners' : ("""vertices in
                            ((x < %.16e) & (y < %.16e) & (z < %.16e))
                          | ((x > %.16e) & (y < %.16e) & (z < %.16e))
                          | ((x > %.16e) & (y > %.16e) & (z < %.16e))
                          | ((x < %.16e) & (y > %.16e) & (z < %.16e))
                          | ((x < %.16e) & (y < %.16e) & (z > %.16e))
                          | ((x > %.16e) & (y < %.16e) & (z > %.16e))
                          | ((x > %.16e) & (y > %.16e) & (z > %.16e))
                          | ((x < %.16e) & (y > %.16e) & (z > %.16e))
                          """ % ( lbnx, lbny, lbnz,
                                  rtfx, lbny, lbnz,
                                  rtfx, rtfy, lbnz,
                                  lbnx, rtfy, lbnz,
                                  lbnx, lbny, rtfz,
                                  rtfx, lbny, rtfz,
                                  rtfx, rtfy, rtfz,
                                  lbnx, rtfy, rtfz ), 'vertex'),
        }
    else:
        lbnx, lbny, = lbn
        rtfx, rtfy, = rtf
        dx = abs(rtfx-lbnx)
        dy = abs(rtfy-lbny)
        lbnx, lbny = (lbnx+dx*eps, lbny+dy*eps,)
        rtfx, rtfy = (rtfx-dx*eps, rtfy-dy*eps,)
        regions = {
            'Bottom' : ('vertices in (y < %.16e)' % lbny, kind),
            'Top' : ('vertices in (y > %.16e)' % rtfy, kind),
            'Left' : ('vertices in (x < %.16e)' % lbnx, kind),
            'Right' : ('vertices in (x > %.16e)' % rtfx, kind),
            'Corners' : ("""vertices in
                              ((x < %.16e) & (y < %.16e))
                            | ((x > %.16e) & (y < %.16e))
                            | ((x > %.16e) & (y > %.16e))
                            | ((x < %.16e) & (y > %.16e))
                            """ % ( lbnx, lbny,
                                    rtfx, lbny,
                                    rtfx, rtfy,
                                    lbnx, rtfy ), 'vertex'),
        }

    return regions

def get_box_volume(dim, lbn, rtf=None):
    """Volume of a box aligned with coordinate axes.

    Parameters:

    dim : int
        Space dimension
    lbn : tuple
        Left bottom near point coordinates if rtf is not None. If rtf is
        None, lbn are the (positive) distances from the origin.
    rtf : tuple
        Right top far point coordinates.

    Returns:

    volume : float
        The box volume.
    """
    if rtf is None:
        lbn, rtf = -nm.array(lbn), lbn

    if dim == 3:
        lbnx, lbny, lbnz = lbn
        rtfx, rtfy, rtfz = rtf
        return abs(rtfx-lbnx)*abs(rtfy-lbny)*abs(rtfz-lbnz)
    else:
        lbnx, lbny, = lbn
        rtfx, rtfy, = rtf
        return abs(rtfx-lbnx)*abs(rtfy-lbny)

def get_lattice_volume(axes):
    r"""
    Volume of a periodic cell in a rectangular 3D (or 2D) lattice.

    Parameters
    ----------
    axes : array
        The array with the periodic cell axes :math:`a_1, \dots, a_3` as rows.

    Returns
    -------
    volume : float
        The periodic cell volume :math:`V = (a_1 \times a_2) \cdot a_3`. In 2D
        :math:`V = |(a_1 \times a_2)|` with zeros as the third components of
        vectors :math:`a_1`, :math:`a_2`.
    """
    axes = nm.asarray(axes)

    dim = axes.shape[0]

    if dim == 2:
        volume = nm.abs(nm.cross(axes[0], axes[1]))

    elif dim == 3:
        volume = nm.dot(nm.cross(axes[0], axes[1]), axes[2])

    else:
        raise ValueError('wrong axes shape! (%s)' % axes.shape)

    return volume

def get_volume(problem, field_name, region_name, quad_order=1):
    """
    Get volume of a given region using integration defined by a given
    field. Both the region and the field have to be defined in
    `problem`.
    """
    from sfepy.discrete import FieldVariable

    field = problem.fields[field_name]
    var = FieldVariable('u', 'parameter', field, 1,
                        primary_var_name='(set-to-None)')

    vol = problem.evaluate('ev_volume.%d.%s( u )' % (quad_order, region_name),
                           u=var)

    return vol

def set_nonlin_states(variables, nl_state, problem):
    """
    Setup reference state for nonlinear homogenization

    Parameters
    ----------
    variables : dict
        All problem variables
    nl_state : reference state
    problem : problem description
    """

    if nl_state is not None:
        var_names = nl_state['variables']
        var_fun = nl_state['set_states']
        pvar_names = []
        for ivar in var_names:
            if ivar in variables:
                pvar_names.append(ivar)
        states = var_fun(problem, pvar_names, variables)
        for ivar in pvar_names:
            variables[ivar].set_data(states[ivar])

def rm_multi(s):
    idx = s.rfind('|multiprocessing_')
    return s[:idx] if idx > 0 else s
