from sfepy.base.base import *
from sfepy.homogenization.utils import coor_to_sym

shared = Struct()

#
# TODO : interpolate fvars to macro times. ?mid-points?
#

def convolve_field_scalar( fvars, pvars, iel, ts ):
    """\int_o^t f(t-s) p(s) ds, t is given by step

    f: fvars ... field variables, defined in a micro domain
    p: pvars ... point variables, a scalar in a point of macro-domain, FMField
    style

    fvars, pvars are instances of History(-duck).
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

def convolve_field_sym_tensor( fvars, pvars, dim, iel, ts ):

    step0 = max( 0, ts.step - fvars[0,0]['u'].steps[-1] )
    print step0, ts.step

    val = nm.zeros_like( fvars[0,0]['u'][0] )
    for ik in xrange( step0, ts.step + 1 ):
        print ' ', ik, ts.step-ik
        for ir in range( dim ):
            for ic in range( dim ):
                ii = coor_to_sym( ir, ic, dim )
                vf = fvars[ir,ic]['u'][ts.step-ik]
                vp = pvars[ik][iel,0,ii,0]
                val += vf * vp * ts.dt
    return val

def compute_u_corr_steady( corrs_rs, strain, corrs_pressure, pressure,
                           dim, iel ):
    """iel = element number"""
    u_corr = corrs_pressure['u'].data * pressure[iel,0,0,0]
    for ir in range( dim ):
        for ic in range( dim ):
            ii = coor_to_sym( ir, ic, dim )
            u_corr += corrs_rs[ir,ic]['u'].data * strain[iel,0,ii,0]
    return u_corr

def compute_u_corr_time( corrs_rs, dstrains, corrs_pressure, pressures,
                         dim, iel, ts ):
    u_corr = convolve_field_scalar( corrs_pressure['u'], pressures,
                                    iel, ts )
    u_corr += convolve_field_sym_tensor( corrs_rs, dstrains,
                                         dim, iel, ts )
    return u_corr


def recover_bones( problem, ts, strain, dstrains, pressure, pressures,
                   corrs_rs, corrs_pressure,
                   corrs_time_rs, corrs_time_pressure ):

    print strain
    print strain.shape
    print dstrains
    print pressure
    print pressure.shape
    print pressures

    dim = problem.domain.mesh.dim
    u_corr_steady = compute_u_corr_steady( corrs_rs, strain,
                                           corrs_pressure, pressure,
                                           dim, 0 )
    print u_corr_steady

    u_corr_time = compute_u_corr_time( corrs_time_rs, dstrains,
                                       corrs_time_pressure, pressures,
                                       dim, 0, ts )
    print u_corr_time
