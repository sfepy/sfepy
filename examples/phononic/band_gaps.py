# c: 25.09.2007, r: 16.12.2008
import os
import numpy as nm

from sfepy import data_dir
from sfepy.fem import MeshIO
from sfepy.homogenization.utils import get_volume

import coef_conf_elastic as cconf
from parametric import vary_incident_wave_dir

filename_mesh = data_dir + '/meshes/2d/special/circle_in_square.mesh'
## filename_mesh = data_dir + '/meshes/2d/special/circle_in_square_small.mesh'
## filename_mesh = data_dir + '/meshes/3d/special/cube_sphere.mesh'
## filename_mesh = data_dir + '/meshes/3d/special/cube_cylinder.mesh'

cwd = os.path.split( os.path.join( os.getcwd(), __file__ ) )[0]

homogeneous = False
fig_suffix = '.pdf'

conf_dir = os.path.dirname(__file__)
io = MeshIO.any_from_filename(filename_mesh,
                              prefix_dir=conf_dir)
bbox, dim = io.read_bounding_box(ret_dim = True)

if homogeneous:
    matrix_region = 'Y'
else:
    matrix_region = 'Y1'

options = {
    'save_eig_vectors' : (10, 0),
    'eig_range' : (0, 30), # -> freq_range = eigs[slice(*eig_range)][[0, -1]]
    'freq_margins' : (10, 10), # % of freq_range

    'feps' : 1e-12, # frequency
    'zeps' : 1e-12, # zero finding
    'teps' : 1e-2, # eigenmomentum threshold
    'teps_rel' : True, # eigenmomentum threshold is relative w.r.t. largest one
    'freq_step' : 0.0001, # % of freq_range
#    'eig_vector_transform' : ('select_in_plane', 'z', 1e-1),
    'plot_transform_angle' : None,
    'plot_transform_wave' : ('clip_sqrt', (0, 30)),
    'plot_transform' : ('normalize', (-2, 2)),

    'output_dir' : os.path.join( cwd, 'output/' ),

    'fig_name' : os.path.join( cwd, 'output', 'band_gaps.pdf' ),
    'fig_name_angle' : os.path.join( cwd, 'output', 'band_gaps_angle.pdf' ),
    'fig_name_wave' : os.path.join( cwd, 'output', 'band_gaps_wave.pdf' ),
    
#    'method' : 'eig.sgscipy',

    'eigenmomentum' : {'var' : 'u',
                       'regions' : ['Y2'],
                       'term' : '%.12e * ev_volume_integrate.2.%s( u )'},

    # Used to compute average density.
    'region_to_material' : {'Y1' : 'matrix',
                            'Y2' : 'inclusion',},
    'tensor_names' : {'elastic' : 'D',},
    'volume' : lambda problem, region_name: get_volume(problem,
                                                       'displacement_Y',
                                                       region_name,
                                                       quad_order=1),
    'eig_problem' : 'simple',

    'dispersion' : 'simple',
    'incident_wave_dir' : [1.0, 1.0],
    'dispersion_conf' : {
        'input' : cconf.define_input(filename_mesh,
                                     matrix_region, bbox),
        'module' : cconf,
    },

    'homogeneous' : homogeneous,
    'fig_suffix' : fig_suffix,
#    'parametric_hook' : 'vary_incident_wave_dir',

    'plot_options' : {
        'show' : True,
        'legend' : True,
    },
    'plot_rsc' : {
        'eig_min' : {'linewidth' : 0.5, 'color' : 'b', 'linestyle' : '--'},
        'eig_max' : {'linewidth' : 0.5, 'color' : 'b', 'linestyle' : '-'},
        'params' : {'axes.labelsize': 'large',
                    'text.fontsize': 'large',
                    'legend.fontsize': 'large',
                    'xtick.labelsize': 'large',
                    'ytick.labelsize': 'large',
                    'text.usetex': False},
    },
}

regions = {
    'Y' : ('all', {}),
    'Y1' : ('elements of group 1', {}),
    'Y2' : ('elements of group 2', {}),
    'Y2_Surface': ('r.Y1 *n r.Y2', {'can_cells' : False}),
}

def get_pars( lam, mu, dim, full = False ):
    from sfepy.mechanics.matcoefs import stiffness_tensor_lame, TransformToPlane
    
    if full:
        c = stiffness_tensor_lame( 3, lam, mu )
        if dim == 2:
            tr = TransformToPlane()
            try:
                c = tr.tensor_plane_stress( c3 = c )
            except:
                sym = (dim + 1) * dim / 2
                c = nm.zeros( (sym, sym), dtype = nm.float64 )
        return c
    else:
        return lam, mu

material_1 = {
    'name' : 'matrix',

    # aluminium, in 1e+10 Pa
    'values' : {
        'lam' : get_pars( 5.898, 2.681, dim )[0],
        'mu' : get_pars( 5.898, 2.681, dim )[1],
        'D' : get_pars( 5.898, 2.681, dim, full = True ),
        'density' : 0.2799, # in 1e4 kg/m3
    },
    'flags' : {'special_constant' : True},
}

material_2 = {
    'name' : 'inclusion',

    # epoxy, in 1e+10 Pa
    'values' : {
        'lam' : get_pars( 0.1798, 0.148, dim )[0],
        'mu' : get_pars( 0.1798, 0.148, dim )[1],
        'D' : get_pars( 0.1798, 0.148, dim, full = True ),
        'density' : 0.1142, # in 1e4 kg/m3
    },
    'flags' : {'special_constant' : True},
}

if homogeneous:
    material_2['values'] = material_1['values']

field_0 = {
    'name' : 'displacement_Y',
    'dtype' : nm.float64,
    'shape' : dim,
    'region' : 'Y',
    'approx_order' : 1,
}

field_1 = {
    'name' : 'displacement_Y2',
    'dtype' : nm.float64,
    'shape' : dim,
    'region' : 'Y2',
    'approx_order' : 1,
}

variables = {
    'u' : ('unknown field', 'displacement_Y2', 0),
    'v' : ('test field', 'displacement_Y2', 'u'),
}

ebc_1 = {
    'name' : 'ZeroSurface',
    'region' : 'Y2_Surface',
    'dofs' : {'u.all' : 0.0},
}

equations = {
    'lhs' : """dw_lin_elastic_iso.2.Y2( inclusion.lam, inclusion.mu, v, u )""",
    'rhs' : """dw_mass_vector.2.Y2( inclusion.density, v, u )""",
}

def clip( data, plot_range ):
    return nm.clip( data, *plot_range )

def clip_sqrt( data, plot_range ):
    return nm.clip( nm.sqrt( data ), *plot_range )

def normalize( data, plot_range ):
    aux = nm.arctan( data )
    return clip( aux, plot_range )

##
# 02.10.2007, c
def select_in_plane( vec, shape, normal_direction, eps ):
    n_nod, dim = shape
    dir_vecs = {2 : {'x': 0, 'y' : 1, 'z' : 1},
               3 : {'x': 0, 'y' : 1, 'z' : 2}}
    ident = nm.eye( dim, dtype = nm.float64 )
    dir_vec = ident[:,dir_vecs[dim][normal_direction]]

    proj = nm.dot( nm.reshape( vec, (n_nod, dim) ), dir_vec )
    if nm.any( nm.abs( proj ) > eps ):
        return nm.zeros_like( vec ), True
    else:
        return vec, False
