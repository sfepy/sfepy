import numpy as nm

import sfepy.discrete.fem.periodic as per
from sfepy.discrete.fem.mesh import Mesh
from sfepy.mechanics.matcoefs import stiffness_from_youngpoisson
from sfepy.homogenization.utils import define_box_regions
import sfepy.homogenization.coefs_base as cb
from sfepy import data_dir

# material function
def get_mat(coors, mode, pb):
    if mode == 'qp':
        cnf = pb.conf
        # get material coefficients
        if hasattr(cnf, 'opt_data'):
            # from optim.
            E_f, nu_f, E_m, nu_m  = cnf.opt_data['mat_params']
        else:
            # given values
            E_f, nu_f, E_m, nu_m  = 160.e9, 0.28, 5.e9, 0.45

        nqp = coors.shape[0]
        nel = pb.domain.mesh.n_el
        nqpe = nqp // nel
        out = nm.zeros((nqp, 6, 6), dtype=nm.float64)

        # set values - matrix
        D_m = stiffness_from_youngpoisson(3, E_m, nu_m)
        Ym = pb.domain.regions['Ym'].get_cells()
        idx0 = (nm.arange(nqpe)[:,nm.newaxis] * nm.ones((1, Ym.shape[0]),
                    dtype=nm.int32)).T.flatten()
        idxs = (Ym[:,nm.newaxis] * nm.ones((1, nqpe),
                    dtype=nm.int32)).flatten() * nqpe
        out[idxs + idx0,...] = D_m

        # set values - fiber
        D_f = stiffness_from_youngpoisson(3, E_f, nu_f)
        Yf = pb.domain.regions['Yf'].get_cells()
        idx0 = (nm.arange(nqpe)[:,nm.newaxis] * nm.ones((1, Yf.shape[0]),
                    dtype=nm.int32)).T.flatten()
        idxs = (Yf[:,nm.newaxis] * nm.ones((1, nqpe),
                    dtype=nm.int32)).flatten() * nqpe
        out[idxs + idx0,...] = D_f

        return {'D': out}

def optimization_hook(pb):
    cnf = pb.conf
    out = []
    yield pb, out

    if hasattr(cnf, 'opt_data'):
        # store homogenized tensor
        pb.conf.opt_data['D_homog'] = out[-1].D.copy()

    yield None

def define(is_opt=False):
    filename_mesh = data_dir + '/meshes/3d/matrix_fiber_rand.vtk'

    mesh = Mesh.from_file(filename_mesh)
    bbox = mesh.get_bounding_box()

    regions = {
        'Y' : 'all',
        'Ym' : ('cells of group 7', 'cell'),
        'Yf' : ('r.Y -c r.Ym', 'cell'),
    }

    regions.update(define_box_regions(3, bbox[0], bbox[1]))

    functions = {
        'get_mat': (lambda ts, coors, mode=None, problem=None, **kwargs:
                    get_mat(coors, mode, problem),),
        'match_x_plane' : (per.match_x_plane,),
        'match_y_plane' : (per.match_y_plane,),
        'match_z_plane' : (per.match_z_plane,),
    }

    materials = {
        'mat': 'get_mat',
    }

    fields = {
        'corrector' : ('real', 3, 'Y', 1),
    }

    variables = {
        'u': ('unknown field', 'corrector'),
        'v': ('test field', 'corrector', 'u'),
        'Pi': ('parameter field', 'corrector', 'u'),
        'Pi1': ('parameter field', 'corrector', '(set-to-None)'),
        'Pi2': ('parameter field', 'corrector', '(set-to-None)'),
    }

    ebcs = {
        'fixed_u' : ('Corners', {'u.all' : 0.0}),
    }

    epbcs = {
        'periodic_x' : (['Left', 'Right'], {'u.all' : 'u.all'}, 'match_x_plane'),
        'periodic_y' : (['Near', 'Far'], {'u.all' : 'u.all'}, 'match_y_plane'),
        'periodic_z' : (['Top', 'Bottom'], {'u.all' : 'u.all'}, 'match_z_plane'),
    }

    all_periodic = ['periodic_%s' % ii for ii in ['x', 'y', 'z'][:3]]

    options = {
        'coefs': 'coefs',
        'requirements': 'requirements',
        'volume': { 'variables' : ['u'], 'expression' : 'ev_volume.5.Y( u )' },
        'output_dir': 'output',
        'coefs_filename': 'coefs_le',
    }

    equation_corrs = {
        'balance_of_forces':
            """dw_lin_elastic.5.Y(mat.D, v, u)
                = - dw_lin_elastic.5.Y(mat.D, v, Pi)"""
    }

    coefs = {
        'D' : {
            'requires' : ['pis', 'corrs_rs'],
            'expression' : 'dw_lin_elastic.5.Y(mat.D, Pi1, Pi2 )',
            'set_variables': [('Pi1', ('pis', 'corrs_rs'), 'u'),
                              ('Pi2', ('pis', 'corrs_rs'), 'u')],
            'class' : cb.CoefSymSym,
        },
        'vol': {
            'regions': ['Ym', 'Yf'],
            'expression': 'ev_volume.5.%s(u)',
            'class': cb.VolumeFractions,
            },
        'filenames' : {},
    }

    requirements = {
        'pis' : {
            'variables' : ['u'],
            'class' : cb.ShapeDimDim,
        },
        'corrs_rs' : {
            'requires' : ['pis'],
            'ebcs' : ['fixed_u'],
            'epbcs' : all_periodic,
            'equations' : equation_corrs,
            'set_variables' : [('Pi', 'pis', 'u')],
            'class' : cb.CorrDimDim,
            'save_name' : 'corrs_le',
            'is_linear': True,
        },
    }

    solvers = {
        'ls' : ('ls.auto_direct', {'use_presolve' : True}),
        'newton' : ('nls.newton', {
            'i_max' : 1,
            'eps_a' : 1e-4,
            'problem': 'linear',
        })
    }

    if is_opt:
        options.update({
            'parametric_hook': 'optimization_hook',
            'float_format': '%.16e',
        })

    return locals()
