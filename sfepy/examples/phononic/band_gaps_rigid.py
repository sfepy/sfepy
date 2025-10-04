"""
Acoustic band gaps in a strongly heterogeneous elastic body with a rigid
inclusion, detected using homogenization techniques.

A reference periodic cell contains three domains: the stiff matrix :math:`Y_m`
and the soft inclusion :math:`Y_c` enclosing the rigid heavy sub-inclusion
:math:`Y_r`.
"""
import numpy as nm

from sfepy import data_dir
from sfepy.base.base import Struct
from sfepy.base.ioutils import InDir
from sfepy.discrete.fem import extend_cell_data
from sfepy.linalg import norm_l2_along_axis
from sfepy.homogenization.coefficients import Coefficients

from sfepy.examples.phononic.band_gaps_conf import (BandGapsRigidConf,
                                                    get_pars, clip)

clip # Make pyflakes happy...

incwd = InDir(__file__)

dim = 2

if dim == 3:
    filename = data_dir + '/meshes/3d/special/cube_sphere.mesh'

else:
    filename = data_dir + '/meshes/2d/special/circle_in_square.mesh'


output_dir = incwd('output/band_gaps_rigid')

# Rigid inclusion diameter.
yr_diameter = 0.125

# aluminium, SI units
D_m = get_pars(2, 5.898e10, 2.681e10)
density_m = 2799.0

# epoxy, SI units
D_c = get_pars(2, 1.798e9, 1.48e9)
density_c = 1142.0

# lead, SI units, does not matter
D_r = get_pars(dim, 4.074e10, 5.556e9)
density_r = 11340.0

mat_pars = Coefficients(D_m=D_m, density_m=density_m,
                        D_c=D_c, density_c=density_c,
                        D_r=D_r, density_r=density_r)

region_selects = Struct(matrix='cells of group 1',
                        inclusion='cells of group 2')

corrs_save_names = {'evp' : 'evp'}

evp_options = {
    'eigensolver' : 'eig.sgscipy',
    'save_eig_vectors' : (12, 0),
    'scale_epsilon' : 1.0,
    'elasticity_contrast' : 1.0,
}

eigenmomenta_options = {
    # eigenmomentum threshold,
    'threshold' : 1e-1,
    # eigenmomentum threshold is relative w.r.t. largest one,
    'threshold_is_relative' : True,
}

band_gaps_options = {
    'fixed_freq_range' : (0., 35000.), # overrides eig_range!

    'freq_eps' : 1e-7, # frequency
    'zero_eps' : 1e-12, # zero finding
    'freq_step' : 0.01, # % of freq_range

    'log_save_name' : 'band_gaps.log',
    'raw_log_save_name' : 'raw_eigensolution.npz',
}

options = {
    'post_process_hook' : 'post_process',

    'plot_transform' : ('clip', (-7000, 7000)),

    'fig_name' : 'band_gaps',
    'fig_suffix' : '.pdf',

    'coefs_filename' : 'coefs.txt',

    'plot_options' : {
        'show' : True, # Show figure.
        'legend' : True, # Show legend.
    },
    'float_format' : '%.16e',
}

def select_yr_circ(coors, diameter=None):
    r = norm_l2_along_axis(coors)
    out = nm.where(r < diameter)[0]

    if out.shape[0] <= 3:
        raise ValueError('too few nodes selected! (%d)' % out.shape[0])

    return out

def _select_yr_circ(coors, domain=None, diameter=None):
    return select_yr_circ(coors, diameter=yr_diameter)

def post_process(out, problem, mtx_phi):
    for key in list(out.keys()):
        ii = int(key[1:])
        vec = mtx_phi[:,ii].copy()
        problem.set_default_state(vec)

        strain = problem.evaluate('ev_cauchy_strain.i.Y_c(u)',
                                  verbose=False, mode='el_avg')
        strain = extend_cell_data(strain, problem.domain, 'Y_c')
        out['strain%03d' % ii] = Struct(name='output_data',
                                        mode='cell', data=strain,
                                        dofs=None)
    return out

conf = BandGapsRigidConf(filename, 1, region_selects, mat_pars, options,
                         evp_options, eigenmomenta_options, band_gaps_options,
                         corrs_save_names=corrs_save_names, incwd=incwd,
                         output_dir=output_dir, select_yr=_select_yr_circ)

define = lambda: conf.conf.to_dict()
