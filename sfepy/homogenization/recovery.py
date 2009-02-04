from sfepy.base.base import *
from sfepy.homogenization.utils import coor_to_sym

shared = Struct()

#
# TODO : interpolate fvars to macro times. ?mid-points?
#

def convolve_field_scalar( fvars, pvars, iel, ts ):
    """
    \int_0^t f(t-s) p(s) ds, t is given by step

    f: fvars ... scalar field variables, defined in a micro domain
    p: pvars ... scalar point variables, a scalar in a point of macro-domain,
    FMField style

    pvars have shape [step][fmf dims]
    fvars have shape [n_step][var dims]
    """

    step0 = max( 0, ts.step - fvars.steps[-1] )
    print step0, ts.step

    val = nm.zeros_like( fvars[0] )
    for ik in xrange( step0, ts.step + 1 ):
        print ' ', ik, ts.step-ik
        vf = fvars[ts.step-ik]
        vp = pvars[ik][iel,0,0,0]
        val += vf * vp * ts.dt

    return val

def convolve_field_sym_tensor( fvars, pvars, var_name, dim, iel, ts ):
    """
    \int_0^t f^{ij}(t-s) p_{ij}(s) ds, t is given by step

    f: fvars ... field variables, defined in a micro domain
    p: pvars ... sym. tensor point variables, a scalar in a point of
    macro-domain, FMField style

    pvars have shape [step][fmf dims]
    fvars have shape [dim, dim][var_name][n_step][var dims]
    """

    step0 = max( 0, ts.step - fvars[0,0][var_name].steps[-1] )
    print step0, ts.step

    val = nm.zeros_like( fvars[0,0][var_name][0] )
    for ik in xrange( step0, ts.step + 1 ):
        print ' ', ik, ts.step-ik
        for ir in range( dim ):
            for ic in range( dim ):
                ii = coor_to_sym( ir, ic, dim )
                vf = fvars[ir,ic][var_name][ts.step-ik]
                vp = pvars[ik][iel,0,ii,0]
                val += vf * vp * ts.dt
    return val

def compute_u_corr_steady( corrs_rs, strain, corrs_pressure, pressure,
                           dim, iel ):
    """
    $\sum_{ij}\left [ \bar\omegabf^{ij}e_{ij}(\ub(t)) \right ]
    + \omegabf^{*,P} p(t)$

    iel = element number
    """
    u_corr = corrs_pressure['u'].data * pressure[iel,0,0,0]
    for ir in range( dim ):
        for ic in range( dim ):
            ii = coor_to_sym( ir, ic, dim )
            u_corr += corrs_rs[ir,ic]['u'].data * strain[iel,0,ii,0]
    return u_corr

def compute_u_corr_time( corrs_rs, dstrains, corrs_pressure, pressures,
                         dim, iel, ts ):
    """
    $\sum_{ij}\left [ \int_0^t \tilde\omegabf^{ij}(t-s)
    \dt{}{s}e_{ij}(\ub(s))\,ds\right ]
    + \int_0^t \tilde\omegabf^P(t-s)\,p(s)\,ds\right ]$
    """
    u_corr = convolve_field_scalar( corrs_pressure['u'], pressures,
                                    iel, ts )
    u_corr += convolve_field_sym_tensor( corrs_rs, dstrains, 'u',
                                         dim, iel, ts )
    return u_corr

def compute_p_corr_steady( corrs_pressure, pressure, iel ):
    """
    $\tilde\pi^P(0)p(t)$
    """
    p_corr = corrs_pressure['p'].data * pressure[iel,0,0,0]
    return p_corr

def compute_p_corr_time( corrs_rs, dstrains, corrs_pressure, pressures,
                         dim, iel, ts ):
    """
    $\sum_{ij} \int_0^t \dt{}{t}\tilde\pi^{ij}(t-s) \dt{}{s}e_{ij}(\ub(s))\,ds
    + \int_0^t \dt{}{t}\tilde\pi^P(t-s)\,p(s)\,ds$
    """
    p_corr = convolve_field_scalar( corrs_pressure['dp'], pressures,
                                    iel, ts )
    p_corr += convolve_field_sym_tensor( corrs_rs, dstrains, 'dp',
                                         dim, iel, ts )
    return p_corr

def compute_u_from_macro( strain, coor, iel ):
    """
    Macro-induced displacements.
    
    e_{ij}^x(\ub(t))\,y_j
    """
    n_nod, dim = coor.shape
    um = nm.empty( (n_nod * dim,), dtype = nm.float64 )
    for ir in range( dim ):
        for ic in range( dim ):
            ii = coor_to_sym( ir, ic, dim )
            um[ir::dim] = strain[iel,0,ii,0] * coor[:,ic]
    return um

def recover_bones( problem, micro_problem,
                   ts, strain, dstrains, pressure, pressures,
                   corrs_rs, corrs_pressure,
                   corrs_time_rs, corrs_time_pressure ):
    """
    note that \tilde{\pi}^P(0) is in corrs_pressure
    -> from time correctors only 'u', 'dp' are needed.
    """
    
    print strain
    print strain.shape
    print dstrains
    print pressure
    print pressure.shape
    print pressures

    dim = problem.domain.mesh.dim
    var = micro_problem.variables['uc']
    micro_coor = var.field.get_coor()[:,:-1]

    for iel in xrange( len( pressures ) ):
        print 'iel:', iel
        u_corr_steady = compute_u_corr_steady( corrs_rs, strain,
                                               corrs_pressure, pressure,
                                               dim, iel )
        u_corr_time = compute_u_corr_time( corrs_time_rs, dstrains,
                                           corrs_time_pressure, pressures,
                                           dim, iel, ts )

        p_corr_steady = compute_p_corr_steady( corrs_pressure, pressure, iel )

        p_corr_time = compute_p_corr_time( corrs_time_rs, dstrains,
                                           corrs_time_pressure, pressures,
                                           dim, iel, ts )
##     print u_corr_steady
##     print u_corr_time
##     print p_coor_steady
##     print p_corr_time

        u_corr = u_corr_steady + u_corr_time
        p_corr = p_corr_steady + p_corr_time

        u_mic = compute_u_from_macro( strain, micro_coor, iel ) + u_corr

        debug()
    
