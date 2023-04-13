r"""
Piezo-elasticity problem - homogenization of a piezoelectric linear elastic
matrix with embedded metalic electrodes, see [1] for details.

[1] E.Rohan, V.Lukes: Homogenization of the fluid-saturated piezoelectric
porous media. International Journal of Solids and Structures 147, 2018,
pages 110-125. https://doi.org/10.1016/j.ijsolstr.2018.05.017
"""
import os
import numpy as nm

from sfepy import data_dir
from sfepy.mechanics.matcoefs import stiffness_from_youngpoisson
from sfepy.homogenization.utils import coor_to_sym, define_box_regions
from sfepy.discrete.fem.mesh import Mesh
from sfepy.base.base import Struct
import sfepy.discrete.fem.periodic as per
import sfepy.homogenization.coefs_base as cb


# Recover fields at the microscpopic level
def recovery_micro(pb, corrs, macro):
    eps0 = macro['eps0']
    mesh = pb.domain.mesh
    regions = pb.domain.regions
    dim = mesh.dim

    if len(macro['u'].shape) == 2:
        mac_strain_Ymc = mac_strain_Ym = macro['strain'][None, ...]
        mac_u_Ymc = macro['u'][None, ...]
    else:
        Ymc_map = regions['Ymc'].get_entities(0)
        Ym_map = regions['Ym'].get_entities(0)
        mac_strain_Ymc = macro['strain'][Ymc_map, ...]
        mac_strain_Ym = macro['strain'][Ym_map, ...]
        mac_u_Ymc = macro['u'][Ymc_map, ...]

    # deformation
    u1, phi = 0, 0

    for ii in range(2):
        u1 += corrs['corrs_k%d' % ii]['u'] * macro['_phi'][ii]
        phi += corrs['corrs_k%d' % ii]['r'] * macro['_phi'][ii]

    for ii in range(dim):
        for jj in range(dim):
            kk = coor_to_sym(ii, jj, dim)
            sidx = f'_{ii}{jj}'
            phi += corrs['corrs_rs']['r' + sidx] * mac_strain_Ym[:, kk]
            u1 += corrs['corrs_rs']['u' + sidx] * mac_strain_Ymc[:, kk]

    u = mac_u_Ymc[..., 0] + eps0 * u1
    mvar = pb.create_variables(['u', 'r', 'svar'])
    e_mac_Ymc = [None] * macro['strain'].shape[1]

    if mac_strain_Ymc.shape[0] > 1:
        for ii in range(dim):
            for jj in range(dim):
                kk = coor_to_sym(ii, jj, dim)
                mvar['svar'].set_data(mac_strain_Ymc[:, kk])
                aux = pb.evaluate('ev_integrate.i2.Ymc(svar)',
                                  mode='el_avg',
                                  var_dict={'svar': mvar['svar']})
                e_mac_Ymc[kk] = aux.squeeze()

        e_mac_Ymc = nm.vstack(e_mac_Ymc).T[:, nm.newaxis, :, nm.newaxis]
    else:
        e_mac_Ymc = mac_strain_Ymc

    mvar['r'].set_data(phi)
    E_mic = pb.evaluate('ev_grad.i2.Ym(r)',
                        mode='el_avg',
                        var_dict={'r': mvar['r']}) / eps0

    mvar['u'].set_data(u1)
    e_mic = pb.evaluate('ev_cauchy_strain.i2.Ymc(u)',
                        mode='el_avg',
                        var_dict={'u': mvar['u']})
    e_mic += e_mac_Ymc

    out = {
        'u': (u, 'u', 'p'),
        'u1': (u1, 'u', 'p'),
        'e_mic': (e_mic, 'u', 'c'),
        'phi': (phi, 'r', 'p'),
        'E_mic': (E_mic, 'r', 'c'),
    }

    out_struct = {}
    for k, v in out.items():
        out_struct[k] = Struct(name='output_data',
                               mode='cell' if v[2] == 'c' else 'vertex',
                               data=v[0],
                               var_name=v[1],
                               dofs=None)

    return out_struct


# Define the local problems and the homogenized coefficients,
# eps0 is the real size of the reference cell.
def define(eps0=1e-3, filename_mesh='meshes/3d/piezo_mesh_micro.vtk'):
    filename_mesh = os.path.join(data_dir, filename_mesh)
    mesh = Mesh.from_file(filename_mesh)
    bbox = mesh.get_bounding_box()
    regions = define_box_regions(mesh.dim, bbox[0], bbox[1], eps=1e-3)

    regions.update({
        'Ymc': 'all',
        # matrix
        'Ym': 'cells of group 1',
        'Ym_left': ('r.Ym *v r.Left', 'vertex'),
        'Ym_right': ('r.Ym *v r.Right', 'vertex'),
        'Ym_bottom': ('r.Ym *v r.Bottom', 'vertex'),
        'Ym_top': ('r.Ym *v r.Top', 'vertex'),
        'Ym_far': ('r.Ym *v r.Far', 'vertex'),
        'Ym_near': ('r.Ym *v r.Near', 'vertex'),
        'Gamma_ms': ('r.Ym *v r.Yc', 'facet', 'Ym'),
        # conductors
        'Yc': ('r.Yc1 +c r.Yc2', 'cell'),
        'Yc1': 'cells of group 2',
        'Yc2': 'cells of group 3',
        'Gamma_s1': ('r.Ym *v r.Yc1', 'facet', 'Ym'),
        'Gamma_s2': ('r.Ym *v r.Yc2', 'facet', 'Ym'),
    })

    options = {
        'coefs_filename': 'coefs_piezo',
        'volume': {'value': nm.prod(bbox[1] - bbox[0])},
        'coefs': 'coefs',
        'requirements': 'requirements',
        'output_dir': 'output',
        'split_results_by': 'region',
        'absolute_mesh_path': True,
        'multiprocessing': False,
        'recovery_hook': recovery_micro,
    }

    fields = {
        'displacement': ('real', 'vector', 'Ymc', 1),
        'potential': ('real', 'scalar', 'Ym', 1),
        'sfield': ('real', 'scalar', 'Ymc', 1),
    }

    variables = {
        # displacement
        'u': ('unknown field', 'displacement'),
        'v': ('test field', 'displacement', 'u'),
        'Pi_u': ('parameter field', 'displacement', 'u'),
        'U1': ('parameter field', 'displacement', '(set-to-None)'),
        'U2': ('parameter field', 'displacement', '(set-to-None)'),
        # potential
        'r': ('unknown field', 'potential'),
        's': ('test field', 'potential', 'r'),
        'Pi_r': ('parameter field', 'potential', 'r'),
        'R1': ('parameter field', 'potential', '(set-to-None)'),
        'R2': ('parameter field', 'potential', '(set-to-None)'),
        # auxiliary
        'svar': ('parameter field', 'sfield', '(set-to-None)'),
    }

    epbcs = {
        'p_ux': (['Left', 'Right'], {'u.all': 'u.all'}, 'match_x_plane'),
        'p_uy': (['Near', 'Far'], {'u.all': 'u.all'}, 'match_y_plane'),
        'p_uz': (['Bottom', 'Top'], {'u.all': 'u.all'}, 'match_z_plane'),
        'p_rx': (['Ym_left', 'Ym_right'], {'r.0': 'r.0'}, 'match_x_plane'),
        'p_ry': (['Ym_near', 'Ym_far'], {'r.0': 'r.0'}, 'match_y_plane'),
        'p_rz': (['Ym_bottom', 'Ym_top'], {'r.0': 'r.0'}, 'match_z_plane'),
    }

    periodic = {
        'per_u': ['per_u_x', 'per_u_y', 'per_u_z'],
        'per_r': ['per_r_x', 'per_r_y', 'per_r_z'],
    }

    # rescale piezoelectric material parameters
    mat_g_sc, mat_d_sc = (eps0, eps0**2)

    materials = {
        'elastic': ({
            'D': {
                'Ym': nm.array([[1.504, 0.656, 0.659, 0, 0, 0],
                                [0.656, 1.504, 0.659, 0, 0, 0],
                                [0.659, 0.659, 1.455, 0, 0, 0],
                                [0, 0, 0, 0.424, 0, 0],
                                [0, 0, 0, 0, 0.439, 0],
                                [0, 0, 0, 0, 0, 0.439]]) * 1e11,
                'Yc': stiffness_from_youngpoisson(3, 200e9, 0.25)}},),
        'piezo': ({
            'g': nm.array([[0, 0, 0, 0, 11.404, 0],
                           [0, 0, 0, 0, 0, 11.404],
                           [-4.322, -4.322, 17.360, 0, 0, 0]]) / mat_g_sc,
            'd': nm.array([[1.284, 0, 0],
                           [0, 1.284, 0],
                           [0, 0, 1.505]]) * 1e-8 / mat_d_sc},),
    }

    functions = {
        'match_x_plane': (per.match_x_plane,),
        'match_y_plane': (per.match_y_plane,),
        'match_z_plane': (per.match_z_plane,),
    }

    ebcs = {
        'fixed_u': ('Corners', {'u.all': 0.0}),
        'fixed_r': ('Gamma_ms', {'r.all': 0.0}),
        'fixed_r1_s1': ('Gamma_s1', {'r.0': 1.0}),
        'fixed_r0_s1': ('Gamma_s1', {'r.0': 0.0}),
        'fixed_r1_s2': ('Gamma_s2', {'r.0': 1.0}),
        'fixed_r0_s2': ('Gamma_s2', {'r.0': 0.0}),
    }

    integrals = {
        'i2': 2,
    }

    solvers = {
        'ls_d': ('ls.auto_direct', {'use_presolve' : True}),
        'ls_i': ('ls.scipy_iterative', {}),
        'ns_ea6': ('nls.newton', {'eps_a': 1e6, 'eps_r': 1e-3,}),
        'ns_ea0': ('nls.newton', {'eps_a': 1e0, 'eps_r': 1e-3,}),
    }

    coefs = {
        'A1': {
            'status': 'auxiliary',
            'requires': ['pis_u', 'corrs_rs'],
            'expression': 'dw_lin_elastic.i2.Ymc(elastic.D, U1, U2)',
            'set_variables': [('U1', ('corrs_rs', 'pis_u'), 'u'),
                              ('U2', ('corrs_rs', 'pis_u'), 'u')],
            'class': cb.CoefSymSym,
        },
        'A2': {
            'status': 'auxiliary',
            'requires': ['corrs_rs'],
            'expression': 'dw_diffusion.i2.Ym(piezo.d, R1, R2)',
            'set_variables': [('R1', 'corrs_rs', 'r'),
                              ('R2', 'corrs_rs', 'r')],
            'class': cb.CoefSymSym,
        },
        'A': {
            'requires': ['c.A1', 'c.A2'],
            'expression': 'c.A1 + c.A2',
            'class': cb.CoefEval,
        },
        'vol': {
            'regions': ['Ym', 'Yc1', 'Yc2'],
            'expression': 'ev_volume.i2.%s(svar)',
            'class': cb.VolumeFractions,
        },
        'eps0': {
            'requires': [],
            'expression': '%e' % eps0,
            'class': cb.CoefEval,
        },
        'filenames': {},
    }

    requirements = {
        'pis_u': {
            'variables': ['u'],
            'class': cb.ShapeDimDim,
        },
        'pis_r': {
            'variables': ['r'],
            'class': cb.ShapeDim,
        },
        'corrs_rs': {
            'requires': ['pis_u'],
            'ebcs': ['fixed_u', 'fixed_r'],
            'epbcs': ['p_ux', 'p_uy', 'p_uz', 'p_rx', 'p_ry', 'p_rz'],
            'equations': {
                'eq1':
                    """dw_lin_elastic.i2.Ymc(elastic.D, v, u)
                     - dw_piezo_coupling.i2.Ym(piezo.g, v, r)
                   = - dw_lin_elastic.i2.Ymc(elastic.D, v, Pi_u)""",
                'eq2':
                    """
                     - dw_piezo_coupling.i2.Ym(piezo.g, u, s)
                     - dw_diffusion.i2.Ym(piezo.d, s, r)
                     = dw_piezo_coupling.i2.Ym(piezo.g, Pi_u, s)""",
            },
            'set_variables': [('Pi_u', 'pis_u', 'u')],
            'class': cb.CorrDimDim,
            'save_name': 'corrs_rs',
            'solvers': {'ls': 'ls_d', 'nls': 'ns_ea6'},
            'is_linear' : True,
        },
    }

    # define requirements and coefficients related to conductors
    bc_conductors = [
        ['fixed_r1_s1', 'fixed_r0_s2'],  # phi = 1 on S1, phi = 0 on S2
        ['fixed_r1_s2', 'fixed_r0_s1'],  # phi = 0 on S1, phi = 1 on S2
    ]

    for k in range(2):
        sk = '%d' % k

        requirements.update({
            'corrs_k' + sk: {
                'requires': ['pis_r'],
                'ebcs': ['fixed_u'] + bc_conductors[k],
                'epbcs': ['p_ux', 'p_uy', 'p_uz', 'p_rx', 'p_ry', 'p_rz'],
                'equations': {
                    'eq1':
                        """dw_lin_elastic.i2.Ymc(elastic.D, v, u)
                         - dw_piezo_coupling.i2.Ym(piezo.g, v, r)
                         = 0""",
                    'eq2':
                        """
                         - dw_piezo_coupling.i2.Ym(piezo.g, u, s)
                         - dw_diffusion.i2.Ym(piezo.d, s, r)
                         = 0"""
                },
                'class': cb.CorrOne,
                'save_name': 'corrs_k' + sk,
                'solvers': {'ls': 'ls_d', 'nls': 'ns_ea0'},
                'is_linear' : True,
            },
        })

        coefs.update({
            'V1_' + sk: {
                'status': 'auxiliary',
                'requires': ['pis_u', 'corrs_k' + sk],
                'expression': 'dw_lin_elastic.i2.Ymc(elastic.D, U1, U2)',
                'set_variables': [('U1', 'corrs_k' + sk, 'u'),
                                  ('U2', 'pis_u', 'u')],
                'class': cb.CoefSym,
            },
            'V2_' + sk: {
                'status': 'auxiliary',
                'requires': ['pis_u', 'corrs_k' + sk],
                'expression': 'dw_piezo_coupling.i2.Ym(piezo.g, U1, R1)',
                'set_variables': [('R1', 'corrs_k' + sk, 'r'),
                                  ('U1', 'pis_u', 'u')],
                'class': cb.CoefSym,
            },
            'V' + sk: {
                'requires': ['c.V1_' + sk, 'c.V2_' + sk],
                'expression': 'c.V1_%s - c.V2_%s' % (sk, sk),
                'class': cb.CoefEval,
            },
        })

    return locals()
