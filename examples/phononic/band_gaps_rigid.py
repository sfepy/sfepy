# c: 14.12.2007, r: 03.11.2008
import os
import numpy as nm

from sfepy import data_dir
from sfepy.base.base import output, pause, debug, Struct
from sfepy.fem import MeshIO
from sfepy.homogenization.utils import get_volume

from gen_mesh import gen_concentric

is_3D = False
generate_2D = False
fig_suffix = '.pdf'

if is_3D:
    filename_mesh = data_dir + '/meshes/3d/special/cube_sphere.mesh'

    out_groups = [1]
    in_groups = [2]
    diameters_g = None
    tepss_g = nm.logspace( -3, -0.5, 11 )
    default_y3_diameter = 0.1
    diameters_g = nm.linspace( 0.075, 0.26, 11 )
else:
    # This mesh is generated, see below.
    filename_mesh = data_dir + '/meshes/2d/special//circle_in_square.vtk'

    out_groups = [1]
    if generate_2D:
        in_groups, diameters_g = gen_concentric( 'tmp/mesh.geo',
                                                 1., 0.2, 0.08, 0.1, 0.6, 7 )
        diameters_g = nm.array( diameters_g[:-2] ) + 0.001
    else:
        os.system("cp %s tmp/mesh.geo" % filename_mesh.replace('.vtk', '.geo'))
        in_groups = [2]
        diameters_g = None
    tepss_g = nm.logspace( -3, -1, 11 )
    default_y3_diameter = 0.25
    os.system("gmsh -2 tmp/mesh.geo -format mesh")
    os.system("script/mesh_to_vtk.py tmp/mesh.mesh %s" % filename_mesh)
    #pause()

cwd = os.path.split( os.path.join( os.getcwd(), __file__ ) )[0]

options = {
    'save_eig_vectors' : (10, 0),
    # Either:
    'eig_range' : (0, 10), # -> freq_range = eigs[slice(*eig_range)][[0, -1]]
    # Or (this has precedence if not None):
    'fixed_eig_range' : (0., 50.),

    'freq_margins' : (10, 10), # % of freq_range

    'feps' : 1e-4, # frequency
    'zeps' : 1e-10, # zero finding
    'teps' : 1e-1, # eigenmomentum threshold
    'teps_rel' : True, # eigenmomentum threshold is relative w.r.t. largest one
    'freq_step' : 0.02, # in percent of freq_range
#    'eig_vector_transform' : ('select_in_plane', 'z', 1e-1),
#    'plot_transform' : ('clip', (-20, 20)),
    'plot_transform' : ('normalize', (-2, 2)),
#    'plot_transform' : None,

    #############################################
#    'parametric_hook' : 'vary_y3_size',
#    'parametric_hook' : 'vary_teps',
    'post_process_hook' : 'post_process',
    'output_dir' : os.path.join( cwd, 'output/' ),
    #############################################

    'eigenmomentum' : {'var' : 'u',
                       'regions' : ['Y2', 'Y3'],
                       'term' : '%.12e * di_volume_integrate.2.%s( u )'},
    # Used to compute average density.
    'region_to_material' : {'Y1' : 'matrix',
                            'Y2' : 'inclusion',
                            'Y3' : 'rigid',},
    'tensor_names' : {'elastic' : 'lame',},
    'volume' : lambda problem, region_name: get_volume(problem,
                                                       'displacement_Y',
                                                       region_name,
                                                       quad_order=1),
    'eig_problem' : 'simple',

    'fig_name' : os.path.join(cwd, 'output', 'band_gaps_sym_025' + fig_suffix),
    'plot_options' : {
        'show' : True, # Show figure.
        'legend' : False, # Show legend.
    },
    'plot_rsc' : { # Resources for all plots.
        'resonance' : {'linewidth' : 0.5, 'color' : 'k', 'linestyle' : '-' },
        'masked' : {'linewidth' : 0.2, 'color' : 'k', 'linestyle' : ':' },
        'x_axis' : {'linewidth' : 1, 'color' : 'k', 'linestyle' : '-' },
        'eig_min' : {'linewidth' : 1, 'color' : 'k', 'linestyle' : '--' },
        'eig_max' : {'linewidth' : 1, 'color' : 'k', 'linestyle' : '-' },
        'strong_gap' : {'linewidth' : 0, 'facecolor' : (1, 1, 0.5) },
        'weak_gap' : {'linewidth' : 0, 'facecolor' : (1, 1, 1) },
        'propagation' : {'linewidth' : 0, 'facecolor' : (0.5, 1, 0.5) },
##         'strong_gap' : {'linewidth' : 0, 'facecolor' : (0.6, 0.6, 0.6) },
##         'weak_gap' : {'linewidth' : 0, 'facecolor' : (0.8, 0.8, 0.8) },
##         'propagation' : {'linewidth' : 0, 'facecolor' : (1, 1, 1) },
        'params' : {'axes.labelsize': 'x-large',
                    'text.fontsize': 'large',
                    'legend.fontsize': 'large',
                    'xtick.labelsize': 'large',
                    'ytick.labelsize': 'large',
                    'text.usetex': False},
    },
    
}

functions = {
    'select_y3_circ' : (lambda coors, domain=None, diameter=None:
                        select_y3_circ(coors, diameter=default_y3_diameter),),
}

regions = {
    'Y' : ('all', {}),
    'Y1' : (' +e '.join( ('elements of group %d' % ig)
                         for ig in out_groups ), {}),
    'Y23' : (' +e '.join( ('elements of group %d' % ig)
                          for ig in in_groups ), {}),
    'Y3' : ('nodes by select_y3_circ', {}),
    'Y2' : ('r.Y23 -e r.Y3', {}),
    'Y23_Surface': ('r.Y1 *n r.Y23', {'can_cells' : False}),
}

material_1 = {
    'name' : 'matrix',

    # aluminium
    'values' : {
        'lam' : 5.898,
        'mu' : 2.681,
        'density' : 0.2799, # in 1e4 kg/m3
    },
    'flags' : {'special_constant' : True},
}

material_2 = {
    'name' : 'inclusion',

    # epoxy, in 1e+10 Pa
    'values' : {
        'lam' : 0.1798,
        'mu' : 0.148,
        'density' : 0.1142, # in 1e4 kg/m3
    },
    'flags' : {'special_constant' : True},
}

material_3 = {
    'name' : 'rigid',

    # lead, in 1e+10 Pa, does not matter
    'values' : {
        'lam' : 4.074 ,
        'mu' : 0.5556,
        'density' : 1.1340, # in 1e4 kg/m3
    },
    'flags' : {'special_constant' : True},
}

conf_dir = os.path.dirname(__file__)
dim = MeshIO.any_from_filename(filename_mesh,
                               prefix_dir=conf_dir).read_dimension()

field_0 = {
    'name' : 'displacement_Y',
    'dtype' : nm.float64,
    'shape' : dim,
    'region' : 'Y',
    'approx_order' : 1,
}

field_1 = {
    'name' : 'displacement_Y23',
    'dtype' : nm.float64,
    'shape' : dim,
    'region' : 'Y23',
    'approx_order' : 1,
}

variables = {
    'u' : ('unknown field', 'displacement_Y23', 0),
    'v' : ('test field', 'displacement_Y23', 'u'),
}

ebc_1 = {
    'name' : 'ZeroSurface',
    'region' : 'Y23_Surface',
    'dofs' : {'u.all' : 0.0},
}

lcbc_1 = {
    'name' : 'RigidBody',
    'region' : 'Y3',
    'dofs' : {'u.all' : 'rigid'},
}

##
# Eigenproblem equations.
# dw_lin_elastic_iso.i1.Y3( rigid.lame, v, u ) should have no effect!
equations = {
    'lhs' : """dw_lin_elastic_iso.2.Y2( inclusion.lam, inclusion.mu, v, u )
             + dw_lin_elastic_iso.2.Y3( rigid.lam, rigid.mu, v, u )""",
    'rhs' : """dw_mass_vector.2.Y2( inclusion.density, v, u )
             + dw_mass_vector.2.Y3( rigid.density, v, u )""",
}

def clip( data, plot_range ):
    return nm.clip( data, *plot_range )

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

def select_y3_circ(coors, diameter=None):
    r = coors[:,0]**2 + coors[:,1]**2
    if dim == 3:
        r += coors[:,2]**2
    r = nm.sqrt( r )

    out = nm.where(r < diameter)[0]

    if out.shape[0] <= 3:
        raise ValueError( 'too few nodes selected! (%d)' % n )

    return out

def extend_cell_data( data, pb, rname, val = None ):
    n_el = pb.domain.shape.n_el
    if data.shape[0] == n_el: return data

    if val is None:
        if data.shape[2] > 1: # Vector.
            val = nm.amin( nm.abs( data ) )
        else: # Scalar.
            val = nm.amin( data )

    edata = nm.empty( (n_el,) + data.shape[1:], dtype = nm.float64 )
    edata.fill( val )
    region = pb.domain.regions[rname]
    offs = region.get_cell_offsets()
    eoffs = pb.domain.get_cell_offsets()
##     print offs
##     print eoffs
##     print pb.domain.mat_ids_to_i_gs
##     pause()
    for group in pb.domain.iter_groups():
        ig = group.ig
        ii = eoffs[ig]
        if ig in region.igs:
            n_cell = region.shape[ig].n_cell
            ir = offs[ig]
            edata[ii+region.cells[ig]] = data[ir:ir+n_cell]
    return edata

def post_process( out, problem, mtx_phi ):
    var = problem.get_variables()['u']

    for key in out.keys():
        ii = int( key[1:] )
        vec = mtx_phi[:,ii].copy()
        var.data_from_any(vec)

        strain = problem.evaluate('de_cauchy_strain.2.Y23( u )', u=var,
                                  verbose=False)
        strain = extend_cell_data( strain, problem, 'Y23' )
        out['strain%03d' % ii] = Struct(name='output_data',
                                        mode='cell', data=strain,
                                        dofs=None)
    return out

def save_log( filename, bg, log_item ):
    """Saves band gaps, valid flags, eigenfrequencies."""

    fd = open( filename, 'w' )
    freq_range = bg.freq_range_margins
    fd.write( log_item )
    fd.write( 'squared: %s\n' % False )
    fd.write( 'n_zeroed: %d\n' % bg.n_zeroed )
    fd.write( 'n_eigs: %d\n' % bg.n_eigs )
    fd.write( 'f0 f1 flag_min f_min v_min flag_max f_max v_max'
              ' kind\ndesc\n' )
    format = "%f %f %d %f %f %d %f %f %s\n%s\n"

    n_row = len( freq_range ) - 1
    fd.write( '%d\n' % n_row )
    for ir in xrange( n_row ):
        f0, f1 = freq_range[[ir, ir+1]]
        gmin, gmax = bg.gaps[ir]
        fd.write( format % ((f0, f1) + tuple( gmin ) + tuple( gmax )
                            + bg.kinds[ir]) )

    fd.write( 'valid resonance\n' )
    freq_range = bg.freq_range_initial
    n_row = len( freq_range )
    fd.write( '%d\n' % n_row )
    valid_in_range = bg.valid[bg.eig_range]
    for ir in xrange( n_row ):
        fd.write( '%d %f\n' % (valid_in_range[ir], freq_range[ir] ) )
    fd.close()

def vary_teps( problem ):
    """Vary eigenmomentum threshold."""
    from sfepy.solvers.ts import get_print_info

    output.prefix = 'vary_teps:'

    if tepss_g is None:
        tepss = nm.logspace( -3, -1, 11 )
    else:
        tepss = tepss_g
    ofn_trunk, output_dir = problem.ofn_trunk, problem.output_dir
    join = os.path.join

    n_digit, aux, d_format = get_print_info( len( tepss ) + 1 )
    
    for ii, teps in enumerate( tepss ):
        output( 'iteration %d: teps %.2e' % (ii, teps) )
        opts = problem.conf.options

        opts.teps = teps

        opts.plot_options['show'] = False
        opts.fig_name = join( output_dir,
                              (('band_gaps_%s' % d_format)
                               + '_teps_%3.2e' + fig_suffix) % (ii, teps) )
        problem.ofn_trunk = ofn_trunk + '_' + (d_format % ii)

        out = []
        yield problem, out

        evp, bg = out[-1]

        filename = join( output_dir,
                         ('band_gaps_%s.txt' % d_format) % ii )
        log_item = '$10^q$: %f\n' % teps
        save_log( filename, bg, log_item )

        yield None

def vary_y3_size( problem ):
    """Vary size of Y3 inclusion."""
    from sfepy.fem import ProblemDefinition
    from sfepy.solvers.ts import get_print_info
    
    output.prefix = 'vary_y3_size:'

    y3_diameters = [0.2, 0.25, 0.3, 0.35, 0.4]
    if diameters_g is None:
        y3_diameters = nm.linspace( 0.15, 0.45, 16 )
    else:
        y3_diameters = diameters_g
#    y3_diameters = [0.45]
    ofn_trunk, output_dir = problem.ofn_trunk, problem.output_dir
    join = os.path.join

    conf = problem.conf
    cr = conf.get_raw( 'regions' )
    n_digit, aux, d_format = get_print_info( len( y3_diameters ) + 1 )
    for ii, diameter in enumerate( y3_diameters ):
        output( 'iteration %d: diameter %3.2f' % (ii, diameter) )
        opts = problem.conf.options

        cr['Y3'] = ('nodes by select_y3_circ( x, y, z, %.5f )' % diameter, {})
        conf.edit( 'regions', cr )
        problem = ProblemDefinition.from_conf( conf )

        problem.save_regions( join( output_dir, ('regions_' + d_format) % ii ),
                              ['Y2', 'Y3'] )
        for region in problem.domain.regions:
            if not region.has_cells_if_can():
                raise ValueError( 'region %s has no cells!' % region.name )

        opts.plot_options['show'] = False
        opts.fig_name = join( output_dir,
                              (('band_gaps_%s' % d_format)
                               + '_y3_%03.2f' + fig_suffix) % (ii, diameter) )
        problem.ofn_trunk = ofn_trunk + '_' + (d_format % ii)
        out = []
        yield problem, out

        evp, bg = out[-1]

        filename = join( output_dir,
                         ('band_gaps_%s.txt' % d_format) % ii )
        log_item = '$r(Y_3)$: %f\n' % diameter
        save_log( filename, bg, log_item )

        yield None

