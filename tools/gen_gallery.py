#!/usr/bin/env python
"""
Generate the images and rst files for gallery of SfePy examples.

The following steps need to be made to regenerate the documentation with the
updated example files:

1. remove doc/examples/*::

   $ rm -rf doc/examples/*

2. generate the files:

   $ ./tools/gen_gallery.py

3. regenerate the documentation::

   $ python setup.py htmldocs
"""
import sys
sys.path.append('.')
import os
import tempfile
import glob
import re
import shutil
import subprocess
from itertools import chain
from argparse import ArgumentParser, RawDescriptionHelpFormatter

import numpy as nm
import matplotlib.image as image

import sfepy
from sfepy.base.base import (get_default, ordered_iteritems,
                             import_file, output, Struct)
from sfepy.base.ioutils import (ensure_path, locate_files, remove_files,
                                edit_filename)
from sfepy.scripts.resview import pv_plot, get_camera_position
import pyvista as pv

omits = [
    '__init__.py',
]

omit_images = [
    'band_gaps_conf.py',
    'dg_plot_1D.py',
    'example_dg_common.py',
    'homogenization_opt.py', # Is parameterized.
    'linear_elasticity_opt.py', # Is parameterized.
    'linear_homogenization_postproc.py',
    'linear_homogenization_up.py', # Used in linear_elastic_mM.py.
    'material_opt.py', # Very long calculation.
    'nonlinear_homogenization.py',
    'piezo_elasticity_micro.py', # Used in piezo_elasticity_macro.py.
    'quantum_common.py',
]

omit_dirs = [
    re.compile('.*output.*/').match,
    re.compile('.*__pycache__/').match,
]

custom = {
    'acoustics/vibro_acoustic3d.py': {
        '_Gamma0': {'view_2d': True, 'max_plots': 2},
        '_Omega1': {'camera': [45, 55, 0.8]},
        '_Omega2': {'camera': [45, 55, 0.8]},
    },
    'acoustics/acoustics3d.py': {
        '_Omega_1': {
            'camera': [75, 135, 1.4],
            'grid_vector1': [1.2, 0, 0],
        },
        '_Omega_2': {
            'camera': [75, 135, 1],
            'grid_vector1': [1.2, 0, 0],
        },
    },
    'acoustics/helmholtz_apartment.py': {
        '': {
            'fields': ['imag.E:wimag.E:f10%:p0', 'mat_id:p1'],
            'force_view_3d': True,
            'grid_vector1': [1.1, 0, 0],
            'camera_position': [-9.22684,-8.37688,10.7623,
                                1.51653,0.122742,-0.646791,
                                0.472758,0.432784,0.767593],
            'color_map': 'seismic',
        },
    },
    'diffusion/cube.py': {
        '': {'camera': [225, 55, 0.7]},
    },
    'diffusion/laplace_time_ebcs.py': {
        '': {'camera': [225, 55, 0.7]},
    },
    'diffusion/laplace_fluid_2d.py': {
        '': {'fields': ['phi:p0', 'phi:t50:p0']},
    },
    'diffusion/laplace_1d.py': {
        '': {
            'fields': ['t:wt:p0', 't:p0'],
            'force_view_3d': True,
            'camera_position': [0.243787,-1.51098,0.486608,
                                0.5,0,0,
                                0,0,1],
        },
    },
    'dg/advection_1D.py': {
        '': {
            'fields': ['p_modal0:r:wp_modal0:p0', '1:vw:p0',
                       'p_modal1:r:wp_modal1:p1', '1:vw:p1',
                       'p_modal2:r:wp_modal2:p2', '1:vw:p2'],
            'force_view_3d': True,
            'camera_position': [0.300893,-2.41044,0.408743,
                                0.4,0,0.55728,
                                0,0,1],
            'grid_vector2': [0, 0, 55],
            'max_plots': 1,
            'color_map': 'cool',
        },
    },
    'diffusion/laplace_coupling_lcbcs.py': {
        '': {
            'fields': ['u1:wu1:p0', 'u1:vw:p0',
                       'u2:wu2:p1', 'u2:vw:p1'],
            'force_view_3d': True,
        },
    },
    'diffusion/laplace_iga_interactive.py': {
        'command': 'python3 sfepy/examples/diffusion/laplace_iga_interactive.py -o output',
        'result': 'output/concentric_circles.vtk',
        'dim': 2,
        'sfepy-view-options': {
            '': {
            },
        },
    },
    'diffusion/laplace_refine_interactive.py': {
        'command': 'python3 sfepy/examples/diffusion/laplace_refine_interactive.py output',
        'result': 'output/hanging.vtk',
        'sfepy-view-options': {
            '': {
                'fields': ['u:wu', '1:vw'],
                'camera_position': [-1.27046,-1.02066,2.28476,
                                    -0.0525852,0.116151,0.00474463,
                                    0.169164,0.842974,0.510664],
            },
        },
    },
    'diffusion/laplace_shifted_periodic.py': {
        'command': 'python3 sfepy/examples/diffusion/laplace_shifted_periodic.py',
        'result': 'laplace_shifted_periodic.vtk',
        'sfepy-view-options': {
            '': {
                'fields': ['u:wu:f0.5', '1:vw'],
                'camera_position': [-1.42703,0.70237,1.33953,
                                    -0.0817071,-0.0547046,0.0280939,
                                    0.616241,-0.216582,0.757192],
            },
        },
    },
    'diffusion/poisson_iga.py': {
        '': {
            'fields': ['t:wt:p0', 't:vw:o0.4:p0'],
            'force_view_3d': True,
            'camera': [30, 60, 1.],
        },
    },
    'diffusion/poisson_parallel_interactive.py': {
        'command': 'python3 sfepy/examples/diffusion/poisson_parallel_interactive.py output-parallel -2 --shape=101,101',
        'result': 'output-parallel/sol.h5',
        'sfepy-view-options': {
            '': {
                'fields': ['u:wu', '1:vw'],
                'camera_position': [-1.00607,-1.49437,0.843106,
                                    0.0141441,-0.0501477,-0.0173962,
                                    0.309532,0.315757,0.896932],
            },
        },
    },
    'diffusion/poisson_parametric_study.py': {
        'command': 'sfepy-run sfepy/examples/diffusion/poisson_parametric_study.py',
        'result': 'output/r_omega1/circles_in_square_0_1_2_3_4_5_6.vtk',
        'dim': 2,
        'sfepy-view-options': {
            '': {
            },
        },
    },
    'diffusion/sinbc.py': {
        '_t': {
            'fields': ['t:wt:p0', '1:vw:p0'],
            'force_view_3d': True,
            'camera': [0, -45, 1.],
        },
        '_grad': {
            'fields': ['grad:g:f0.01:p0', '1:vw:p0'],
            'force_view_3d': True,
            'camera': [0, -45, 1.5],
        },
    },
    'diffusion/time_heat_equation_multi_material.py': {
        '': {
            'isosurfaces': 10,
            'outline': True,
            'camera': [-50, -230, 1],
        },
    },
    'diffusion/time_poisson_interactive.py': {
        'command': 'python3 sfepy/examples/diffusion/time_poisson_interactive.py -p',
        'image': 'time_poisson_interactive_probe_04.png',
        'result': 'domain.10.vtk',
        'sfepy-view-options': {
            '': {
            },
        },
    },
    'homogenization/linear_homogenization.py': {
        'command': 'sfepy-run sfepy/examples/homogenization/linear_homogenization.py',
        'result': 'output/corrs_le.vtk',
        'sfepy-view-options': {
            '': {
                'max_plots': 3,
                'camera_position': [-8.57147,2.10122,2.76537,
                                    0.37,2.04405,2.14477,
                                    0,0,1],
                'show_labels': False,
                'show_scalar_bars': False,
            }
        },
    },
    'homogenization/perfusion_micro.py': {
        'command': 'sfepy-run sfepy/examples/homogenization/perfusion_micro.py',
        'result': 'output/corrs_3d_2ch.vtk',
        'sfepy-view-options': {
            '_etaA_YM': {
                'camera_position': [-1.69531,-0.966016,2.3648,
                                    0.480252,0.486045,0.46669,
                                    0.51024,0.292658,0.808707],
            },
            '_etaB_YM': {
                'camera_position': [-1.69531,-0.966016,2.3648,
                                    0.480252,0.486045,0.46669,
                                    0.51024,0.292658,0.808707],
            },
            # '_gamma_A_1_YA': {
            #     'camera_position': [-0.96887,-2.21765,1.25047,
            #                         1,0.4,0.5,
            #                         0.107002,0.198813,0.974179],
            #     'grid_vector1': [1.2,0,0],
            # },
            # '_gamma_A_2_YA': {
            #     'camera_position': [-0.96887,-2.21765,1.25047,
            #                         1,0.4,0.5,
            #                         0.107002,0.198813,0.974179],
            #     'grid_vector1': [1.2,0,0],
            # },
            # '_gamma_A_3_YA': {
            #     'camera_position': [-0.96887,-2.21765,1.25047,
            #                         1,0.4,0.5,
            #                         0.107002,0.198813,0.974179],
            #     'grid_vector1': [1.2,0,0],
            # },
            # '_gamma_B_1_YB': {
            #     'camera_position': [-0.96887,-2.21765,1.25047,
            #                         1,0.4,0.5,
            #                         0.107002,0.198813,0.974179],
            #     'grid_vector1': [1.2,0,0],
            # },
            # '_gamma_B_2_YB': {
            #     'camera_position': [-0.96887,-2.21765,1.25047,
            #                         1,0.4,0.5,
            #                         0.107002,0.198813,0.974179],
            #     'grid_vector1': [1.2,0,0],
            # },
            # '_gamma_B_3_YB': {
            #     'camera_position': [-0.96887,-2.21765,1.25047,
            #                         1,0.4,0.5,
            #                         0.107002,0.198813,0.974179],
            #     'grid_vector1': [1.2,0,0],
            # },
            # '_gamma_m_YM': {
            #     'camera_position': [-1.69531,-0.966016,2.3648,
            #                         0.480252,0.486045,0.46669,
            #                         0.51024,0.292658,0.808707],
            # },
            '_gamma_p_YM': {
                'camera_position': [-1.69531,-0.966016,2.3648,
                                    0.480252,0.486045,0.46669,
                                    0.51024,0.292658,0.808707],
            },
            # '_oneA_YM': {
            #     'camera_position': [-1.69531,-0.966016,2.3648,
            #                         0.480252,0.486045,0.46669,
            #                         0.51024,0.292658,0.808707],
            # },
            # '_oneB_YM': {
            #     'camera_position': [-1.69531,-0.966016,2.3648,
            #                         0.480252,0.486045,0.46669,
            #                         0.51024,0.292658,0.808707],
            # },
            # '_piA_YA': {
            #     'max_plots': 2,
            #     'camera_position': [-2.02509,0.148141,3.9141,
            #                         0.913524,1.35355,0.462869,
            #                         0.714983,0.186315,0.673859],
            #     'grid_vector1': [1.2,0,0],
            # },
            # '_piB_YB': {
            #     'max_plots': 2,
            #     'camera_position': [-2.02509,0.148141,3.9141,
            #                         0.913524,1.35355,0.462869,
            #                         0.714983,0.186315,0.673859],
            #     'grid_vector1': [1.2,0,0],
            # },
        },
    },
    'homogenization/rs_correctors.py': {
        'command': 'python3 sfepy/examples/homogenization/rs_correctors.py -n',
        'result': 'corrs_elastic.vtk',
        'sfepy-view-options': {
            '': {
                'max_plots': 2,
                'camera_position': [1.16983,0.660181,7.63496,
                                    1.16983,0.660181,
                                    0,0,1,0],
            },
        },
    },
    'large_deformation/active_fibres.py': {
        'command_0': 'sfepy-run sfepy/examples/large_deformation/active_fibres.py',
        'command_1': 'sfepy-view output/hsphere8_fibres.vtk -f fdir0:t2000:p0 1:vs:o0.1:p0 --no-step-time --no-scalar-bars --no-axes --camera-position=-0.0337972,-0.0337972,0.0147184,0.0125,0.0125,0.0129782,0,0,1 --off-screen -o output/hsphere8-fdir0.png',
        'image_1': 'output/hsphere8-fdir0.png',
        'command_2': 'sfepy-view output/hsphere8_fibres.vtk -f fdir1:t2000:p0 1:vs:o0.1:p0 --no-step-time --no-scalar-bars --no-axes --camera-position=-0.0337972,-0.0337972,0.0147184,0.0125,0.0125,0.0129782,0,0,1 --off-screen -o output/hsphere8-fdir1.png',
        'image_2': 'output/hsphere8-fdir1.png',
        'result': 'output/hsphere8.h5',
        'result_before_images' : True,
        'sfepy-view-options': {
            '': {
                'fields': ['green_strain:wu:f1:p0', '1:vw:o0.3:p0'],
                'step' : 10,
                'camera_position': [-0.0337972,-0.0337972,0.0147184,
                                    0.0125,0.0125,0.0129782,
                                    0,0,1],
            },
        },
    },
    'large_deformation/compare_elastic_materials.py': {
        'command': 'python3 sfepy/examples/large_deformation/compare_elastic_materials.py -n',
        'image': 'pressure_displacement.png',
        'sfepy-view-options': {
        },
    },
    'large_deformation/gen_yeoh_tl_up_interactive.py': {
        'command': 'python3 sfepy/examples/large_deformation/gen_yeoh_tl_up_interactive.py -pn',
        'image': 'gen_yeoh_tl_up_comparison.png',
        'result': 'domain.10.vtk',
        'sfepy-view-options': {
            '': {
                'fields': ['stress:wu:f1:p0', '1:vw:p0'],
                'camera_position': [-1.89068,-2.47529,0.674151,
                                    0.793688,0.209075,0.574171,
                                    0,0,1],
            }
        },
    },
    'large_deformation/hyperelastic_tl_up_interactive.py': {
        'command': 'python3 sfepy/examples/large_deformation/hyperelastic_tl_up_interactive.py -pn',
        'image': 'hyperelastic_tl_up_comparison.png',
        'result': 'domain.10.vtk',
        'sfepy-view-options': {
            '': {
                'fields': ['stress:wu:f1:p0', '1:vw:p0'],
                'camera_position': [-4.04955,-6.37091,3.78116,
                                    2.06774,-0.866196,1.16309,
                                    0.220544,0.208124,0.952914],
            }
        },
    },
    'linear_elasticity/dispersion_analysis.py': {
        'command_0': 'python3 sfepy/examples/linear_elasticity/dispersion_analysis.py meshes/2d/special/circle_in_square.mesh --log-std-waves --eigs-only --no-show',
        'command_1': 'python3 sfepy/scripts/plot_logs.py output/frequencies.txt -o output/frequencies.png -n',
        'image': 'output/frequencies.png',
        'sfepy-view-options': {
        },
    },
    'linear_elasticity/elastic_contact_planes.py': {
        '': {
            'fields': ['u:wu:p0', '1:vw:p0'],
            'camera': [225, 55, 0.7],
        },
    },
    'linear_elasticity/elastic_contact_sphere.py': {
        '': {
            'fields': ['u:wu:p0', '1:vw:p0'],
            'camera': [225, 55, 0.7],
        },
    },
    'linear_elasticity/elastic_shifted_periodic.py': {
        '': {'fields': ['von_mises_stress:r:wu:p0', '1:vw:p0']},
    },
    'linear_elasticity/elastodynamic.py': {
        '': {
            'fields': ['u:wu:f1e3:p0', '1:vw:p0',
                       'cauchy_strain:p1', 'cauchy_stress:p2'],
        },
    },
    'linear_elasticity/elastodynamic_identification.py': {
        'command_0': 'python3 sfepy/examples/linear_elasticity/elastodynamic_identification.py --output-dir=output/edi',
        'image_0': 'output/edi/res00004.png',
        'command_1': 'python3 sfepy/scripts/plot_logs.py output/edi/pars.txt -o output/edi/pars.png -n',
        'image_1': 'output/edi/pars.png',
        'sfepy-view-options': {
        },
    },
    'linear_elasticity/its2D_4.py': {
        'command_0': 'sfepy-run sfepy/examples/linear_elasticity/its2D_4.py',
        'command_1': 'sfepy-probe sfepy/examples/linear_elasticity/its2D_4.py its2D.h5',
        'image_0': 'its2D_0.png',
        'image_1': 'its2D_1.png',
        'result': 'its2D.h5',
        'dim': 2,
        'sfepy-view-options': {
            '': {
            },
        },
    },
    'linear_elasticity/its2D_5.py': {
        'command': 'sfepy-run sfepy/examples/linear_elasticity/its2D_5.py',
        'image_0': 'its2D_probe_line0.png',
        'image_1': 'its2D_probe_line1.png',
        'result': 'its2D.vtk',
        'dim': 2,
        'sfepy-view-options': {
            '': {
            },
        },
    },
    'linear_elasticity/its2D_interactive.py': {
        'command': 'python3 sfepy/examples/linear_elasticity/its2D_interactive.py -p',
        'image_0': 'its2D_interactive_probe_0.png',
        'image_1': 'its2D_interactive_probe_1.png',
        'result': 'its2D_interactive.vtk',
        'dim': 2,
        'sfepy-view-options': {
            '': {
            },
        },
    },
    'linear_elasticity/linear_elastic_iga.py': {
        '': {
            'fields': ['u:wu:p0', '1:vw:p0'],
            'camera': [-45, 55, 1],
        },
    },
    'linear_elasticity/linear_elastic_interactive.py': {
        'command': 'python3 sfepy/examples/linear_elasticity/linear_elastic_interactive.py',
        'result': 'linear_elasticity.vtk',
        'dim': 2,
        'sfepy-view-options': {
            '': {
                'fields': ['u:wu:p0', '1:vw:p0'],
            },
        },
    },
    'linear_elasticity/linear_elastic_probes.py': {
        'command': 'sfepy-run sfepy/examples/linear_elasticity/linear_elastic_probes.py',
        'image_0': 'cylinder_probe_line.png',
        'image_1': 'cylinder_probe_circle.png',
        'image_2': 'cylinder_probe_ray.png',
        'result': 'cylinder.vtk',
        'sfepy-view-options': {
            '': {
                'fields': ['cauchy_stress:wu:p0', '1:vw:p0'],
                'camera_position': [-0.0496419,-0.132205,0.0339214,
                                    0.0495001,0.00309172,-0.00112878,
                                    0.146854,0.145842,0.978348],
            },
        },
    },
    'linear_elasticity/linear_viscoelastic.py': {
        '': {'camera': [225, 75, 0.88]}
    },
    'linear_elasticity/modal_analysis.py': {
        'command': 'python3 sfepy/examples/linear_elasticity/modal_analysis.py  -s 31,31 -n 3',
        'result': 'eigenshapes.vtk',
        'sfepy-view-options': {
            '': {
                'fields': ['strain000:wu000:f10%:p0', '1:vw:p0',
                           'strain001:wu001:f10%:p1', '1:vw:p1',
                           'strain002:wu002:f10%:p2', '1:vw:p2'],
                'camera_position': [1.1,-0.1,5.41561,
                                    1.1,-0.1,
                                    0,0,1,0],
                'grid_vector1': [1.2,0,0],
            },
        },
    },
    'linear_elasticity/modal_analysis_declarative.py': {
        '': {
            'fields': ['u003:wu003:f30%:p0', '1:vw:p0'],
            'camera_position': [-2.30562,-2.2604,-0.325838,
                                -0.0714771,0.0374911,0.0287214,
                                -0.0245554,-0.129086,0.991329],
        },
    },
    'linear_elasticity/multi_node_lcbcs.py': {
        '': {
            'fields': ['u:wu:e'],
            'force_view_3d': True,
            'camera_position': [0.175,0.125,0.735014,
                                0.175,0.125,0,
                                0,1,0],
        },
    },
    'linear_elasticity/multi_point_constraints.py': {
        '': {
            'fields': ['u:wu:f1:p0', '1:vw:p0', 'u:gu:p0'],
            'force_view_3d': True,
            'camera_position': [0.0565,0.0434999,19.8951,
                                0.0565,0.0434999,0,
                                0,1,0],
        },
    },
    'linear_elasticity/seismic_load.py': {
        '': {
            'fields': ['cauchy_stress:wu:f10:p0', '1:vw:p0'],
        },
    },
    'linear_elasticity/shell10x_cantilever.py': {
        '': {
            'fields': ['u_disp:wu_disp:p0', '1:vw:p0',
                       'u_rot:p1', '1:vw:p1'],
            'camera': [-45, 75, 1],
            'grid_vector1': [1, 0, 0],
        },
    },
    'linear_elasticity/shell10x_cantilever_interactive.py': {
        'command': 'python3 sfepy/examples/linear_elasticity/shell10x_cantilever_interactive.py -t bend --plot --no-show output',
        'image': 'output/shell10x_cantilever_convergence_bent.png',
        'result': 'output/shell10x_cantilever.vtk',
        'sfepy-view-options': {
            '': {
                'fields': ['u_disp:wu_disp:p0', '1:vw:p0',
                           'u_rot:p1', '1:vw:p1'],
                'camera': [-45, 75, 1],
                'grid_vector1': [1, 0, 0],
            },
        },
    },
    'linear_elasticity/truss_bridge.py': {
        '': {
            'view_2d': True,
            'fields': ['u:wu:p0', '1:vw:p0', 'S:e:p1'],
            'grid_vector1': [0, 2, 0],
        },
    },
    'linear_elasticity/truss_bridge3d.py': {
        '_Solid': {
            'fields': ['u_solid:wu_solid:f1e3:p0'],
            'camera_position': [-5.912, -6.64883, 1.80888,
                                5.05199, 2.28013, -1.49468,
                                0, 0, 1],
        },
        '_Struct': {
            'fields': ['u_struct:wu_struct:f1e3:p0'],
            'camera_position': [-5.912, -6.64883, 1.80888,
                                5.05199, 2.28013, -1.49468,
                                0, 0, 1],
        },
    },
    'linear_elasticity/two_bodies_contact.py': {
        '': {
            'fields': ['u:wu:f1:p0', '1:vw:wu:f1:p0'],
            'camera_position': [-1.39408,-2.02778,0.937677,
                                -0.0018284,-0.034985,-0.15843,
                                0.208882,0.355227,0.911143],

        }
    },
    'linear_elasticity/wedge_mesh.py': {
        'command': 'sfepy-run sfepy/examples/linear_elasticity/wedge_mesh.py',
        'result': 'beam_w14.vtk',
        'sfepy-view-options': {
            '': {
                'fields': ['u:wu:e:o0.5'],
                'camera_position': [0.927482,-0.574865,0.307926,
                                    0.372897,0.120369,-0.0347131,
                                    -0.326236,0.19556,0.924838],
            }
        },
    },
    'miscellaneous/live_plot.py': {
        'command_0': 'python3 sfepy/examples/miscellaneous/live_plot.py -o output',
        'command_1': 'python3 sfepy/scripts/plot_logs.py output/live_plot.txt -o output/live_plot.png -n',
        'image_1': 'output/live_plot.png',
        'command_2': 'python3 sfepy/scripts/plot_logs.py output/live_plot2.txt -o output/live_plot2.png -n',
        'image_2': 'output/live_plot2.png',
        'sfepy-view-options': {
        },
    },
    'miscellaneous/refine_evp.py': {
        'command': 'python3 sfepy/examples/miscellaneous/refine_evp.py --max-order=5 --max-refine=2 --fig-suffix=.png --no-show',
        'image': 'output/h-refinement-0-laplace-lagrange-primme-none-a.png',
        'sfepy-view-options': {
        },
    },
    'multi_physics/biot_parallel_interactive.py': {
        'command': 'python3 sfepy/examples/multi_physics/biot_parallel_interactive.py output-parallel',
        'result': 'output-parallel/sol.h5',
        'sfepy-view-options': {
            '': {
                'fields': ['u:t1000:p0', '1:vw:o0.1:p0', 'p:p1'],
                'camera_position': [-0.740809,-2.67483,2.43604,
                                    0.579298,-0.106475,0.0404841,
                                    0.183213,0.618527,0.764106],
                'grid_vector1': [1.2, 0, 0]
            },
        },
    },
    'multi_physics/piezo_elasticity.py': {
        '': {'fields': ['u:g:p0', 'cauchy_strain:p1',
                        'elastic_stress:p2', 'piezo_stress:p3'],
             'grid_vector1': [1.2, 0, 0],
             'max_plots': 4},
    },
    'multi_physics/piezo_elastodynamic.py': {
        '': {
            'fields': ['p:wu:f2000:p0', '1:vw:wu:f2000:p0'],
            'camera_position': [-0.0125588,-0.00559266,0.0117482,
                                0.00438669,0.00487109,0.00135715,
                                0.334487,0.333581,0.881387],
            'color_map': 'inferno',
        }
    },
    'multi_physics/thermal_electric.py': {
        'command': 'python3 sfepy/examples/multi_physics/thermal_electric.py',
        'result': 'sfepy/examples/multi_physics/output/circle_in_square.vtk',
        'dim': 2,
        'sfepy-view-options': {
            '.el': {
            },
            '.10': {
            },
        },
    },
    'multi_physics/thermo_elasticity_ess.py': {
        '': {
            'fields': ['T:wu:f1e3:p0', '1:vw:p0'],
            'camera': [-45, 75, 1],
        },
    },
    'multi_physics/thermo_elasticity.py': {
        '': {'camera': [225, 75, 0.9]}
    },
    'navier_stokes/stokes_slip_bc.py': {
        '': {
            'fields': ['u:g:f.25:p0', 'u:o.4:p0', 'p:p1'],
            'camera': [-45, 55, 1],
            'grid_vector1': [0, 1.2, 0]
        },
    },
    'phononic/band_gaps.py': {
        'command': 'sfepy-run sfepy/examples/phononic/band_gaps.py --phonon-band-gaps --phonon-plot -O plot_options={show=False,legend=True},fig_suffix=\'.png\'',
        'image': 'sfepy/examples/phononic/output/band_gaps/band_gaps.png',
        'result': 'sfepy/examples/phononic/output/band_gaps/evp.vtk',
        'sfepy-view-options': {
            '': {
                'camera_position': [1.8,-1.2,7.28269,
                                    1.8,-1.2,0,
                                    0,1,0],
                'grid_vector1': [1.2, 0, 0],
                'grid_vector2': [0, -1.2, 0],
                'max_plots': 4,
                'show_labels': False,
                'show_scalar_bars': False,
            },
        },
    },
    'phononic/band_gaps_rigid.py': {
        'command': 'sfepy-run sfepy/examples/phononic/band_gaps_rigid.py --phonon-band-gaps --phonon-plot -O plot_options={show=False,legend=True},fig_suffix=\'_rigid.png\'',
        'image': 'sfepy/examples/phononic/output/band_gaps_rigid/band_gaps_rigid.png',
        'result': 'sfepy/examples/phononic/output/band_gaps_rigid/evp.vtk',
        'sfepy-view-options': {
            '': {
                'fields': list(chain(
                    *[[f'u{ii:03d}:vs:o.4:p{ii}', f'u{ii:03d}:g:p{ii}']
                      for ii in range(12)]
                )),
                'camera_position': [1.8,-1.2,7.28269,
                                    1.8,-1.2,0,
                                    0,1,0],
                'grid_vector1': [1.2, 0, 0],
                'grid_vector2': [0, -1.2, 0],
                'max_plots': 4,
                'show_labels': False,
                'show_scalar_bars': False,
            },
        },
    },
    'quantum/boron.py': {
        '': {'fields': ['Psi000:p0', 'Psi001:p1', 'Psi002:p2', 'Psi003:p3'],
             'grid_vector1': [1.2, 0, 0],
             'max_plots': 4},
    },
    'quantum/hydrogen.py': {
        '': {'fields': ['Psi000:p0', 'Psi001:p1', 'Psi002:p2', 'Psi003:p3'],
             'grid_vector1': [1.2, 0, 0],
             'max_plots': 4},
    },
    'quantum/oscillator.py': {
        '': {'fields': ['Psi000:p0', 'Psi001:p1', 'Psi002:p2', 'Psi003:p3'],
             'grid_vector1': [1.2, 0, 0],
             'max_plots': 4},
    },
    'quantum/well.py': {
        '': {'fields': ['Psi000:p0', 'Psi001:p1', 'Psi002:p2', 'Psi003:p3'],
             'grid_vector1': [1.2, 0, 0],
             'max_plots': 4},
    },
}


def resview_plot(filename, filename_out, options):
    pv.set_plot_theme("document")
    plotter = pv.Plotter(off_screen=True)

    plotter = pv_plot([filename], options, plotter=plotter, use_cache=False)
    if options.axes_visibility:
        plotter.add_axes(**dict(options.axes_options))

    if options.view_2d:
        plotter.view_xy()
        plotter.show(screenshot=filename_out, window_size=(800, 600))
    else:
        if options.camera_position is not None:
            cpos = nm.array(options.camera_position)
            cpos = cpos.reshape((3, 3))
        elif options.camera:
            zoom = options.camera[2] if len(options.camera) > 2 else 1.
            cpos = get_camera_position(plotter.bounds,
                                       options.camera[0], options.camera[1],
                                       zoom=zoom)
        else:
            cpos = None

        plotter.show(cpos=cpos, screenshot=filename_out, window_size=(800, 600))

def run_resview_plot(*args):
    """
    A fix for the problem that calling :func:`resview_plot()` directly often
    terminates the program.
    """
    from multiprocessing import Process

    process = Process(target=resview_plot, args=args)
    process.start()
    process.join()

def _omit(filename, omits, omit_dirs):
    omit = False

    base = os.path.basename(filename)

    if base in omits:
        omit = True

    for omit_dir in omit_dirs:
        if omit_dir(filename) is not None:
            omit = True
            break

    return omit


def ebase2fbase(ebase):
    return os.path.splitext(ebase)[0].replace(os.path.sep, '-')

def _make_fig_name(fig_base, fig_filename):
    return fig_base + '-' + os.path.basename(fig_filename)

def _get_image_names(custom_options):
    for key, val in custom_options.items():
        if key.startswith('image'):
                yield val

def _get_result_fig_filenames(fig_base, images_dir, custom_view_options):
    suffixes = sorted(custom_view_options.keys())
    for suffix in suffixes:
        yield os.path.join(images_dir, fig_base + suffix + '.png')

def _get_fig_filenames(ebase, images_dir):
    fig_base = ebase2fbase(ebase)

    result_done = False
    if ebase in custom:
        custom_options = custom.get(ebase)
        result_first = custom_options.get('result_before_images', False)
        if 'sfepy-view-options' in custom_options:
            custom_view_options = custom_options['sfepy-view-options']

            if result_first:
                yield from _get_result_fig_filenames(fig_base, images_dir,
                                                     custom_view_options)
                result_done = True

            for fig_filename in _get_image_names(custom_options):
                yield os.path.join(images_dir,
                                   _make_fig_name(fig_base, fig_filename))

        else:
            custom_view_options = custom_options

        if custom_view_options and not result_done:
            yield from _get_result_fig_filenames(fig_base, images_dir,
                                                 custom_view_options)

    else:
        yield os.path.join(images_dir, fig_base + '.png')


def _get_fig_filename(ebase, images_dir, suffix):
    fig_base = ebase2fbase(ebase)

    return os.path.join(images_dir, fig_base + suffix + '.png')


def _make_sphinx_path(path, relative=False):
    if relative:
        aux = path.replace(sfepy.data_dir, '')
        prefix = ('..' + os.path.sep) * aux.count(os.path.sep)
        sphinx_path = prefix[:-1] + aux

    else:
        sphinx_path = path.replace(sfepy.data_dir, '/..')

    return sphinx_path


def _apply_commands(custom_options, ebase, images_dir):
    for key, val in custom_options.items():
        if key.startswith('command'):
            cmd = custom_options[key]
            subprocess.run(cmd.split()).check_returncode()

        if key.startswith('image'):
            shutil.copy(val,
                        os.path.join(images_dir,
                                     _make_fig_name(ebase2fbase(ebase), val)))


def apply_view_options(views, default):
    out = {}
    for kview, view in views.items():
        ov = default.copy()
        for k, v in view.items():
            if k == 'fields':
                fv = [k.split(':') for k in v]
                fv = [(k[0], ':'.join(k[1:])) for k in fv]
                setattr(ov, k, fv)
            else:
                setattr(ov, k, v)

        out[kview] = ov

    return out


def generate_images(images_dir, examples_dir, pattern='*.py'):
    """
    Generate images from results of running examples found in
    `examples_dir` directory.

    The generated images are stored to `images_dir`,
    """
    from sfepy.applications import solve_pde
    from sfepy.solvers.ts_solvers import StationarySolver

    prefix = output.prefix

    output_dir = tempfile.mkdtemp()
    trunk = os.path.join(output_dir, 'result')
    options = Struct(output_filename_trunk=trunk,
                     output_format='vtk',
                     save_ebc=False,
                     save_ebc_nodes=False,
                     save_regions=False,
                     save_regions_as_groups=False,
                     solve_not=False)
    app_options = dict(
        file_format='vtk',
    )

    view_options = Struct(step=0,
                          fields=[], fields_map=[],
                          outline=False,
                          isosurfaces=0,
                          show_edges=False,
                          warp=None,
                          factor=1.,
                          opacity=1.,
                          color_map='viridis',
                          color_limits=None,
                          axes_options=[],
                          axes_visibility=False,
                          grid_vector1=None,
                          grid_vector2=None,
                          max_plots=3,
                          show_labels=False,
                          label_position=[-1, -1, 0, 0.2],
                          scalar_bar_size=[0.15, 0.06],
                          scalar_bar_position=[0.04, 0.92, 0, -1.5],
                          show_scalar_bars=True,
                          camera=[225, 75, 1],
                          camera_position=None,
                          view_2d=False,
                          force_view_3d=False,
                          show_step_time=False)

    ensure_path(images_dir + os.path.sep)

    for ex_filename in locate_files(pattern, examples_dir):
        if _omit(ex_filename, omits + omit_images, omit_dirs):
            continue

        output.level = 0
        output.prefix = prefix
        ebase = ex_filename.replace(examples_dir, '')[1:]
        output('trying "%s"...' % ebase)

        _ebase = ebase.replace(os.path.sep, '/')
        custom_options = custom.get(_ebase)
        if custom_options and 'sfepy-view-options' in custom_options:
            try:
                _apply_commands(custom_options, ebase, images_dir)

                filename = custom_options.get('result')
                dim = custom_options.get('dim')
                custom_view_options = custom_options['sfepy-view-options']

            except subprocess.CalledProcessError:
                filename = None
                output('***** failed! *****')

        else:
            custom_view_options = custom_options

            try:
                problem, state = solve_pde(ex_filename, options=options,
                                           **app_options)
                try:
                    tsolver = problem.get_solver()

                except ValueError:
                    suffix = None

                else:
                    if isinstance(tsolver, StationarySolver):
                        suffix = None

                    else:
                        suffix = tsolver.ts.suffix % (tsolver.ts.n_step - 1)

                filename = problem.get_output_name(suffix=suffix)
                dim = problem.get_dim()

            except KeyboardInterrupt:
                raise

            except:
                filename = None
                output('***** failed! *****')

        if filename is not None:
            if custom_view_options is not None:
                views = apply_view_options(custom_view_options, view_options)
            else:
                views = {'': view_options.copy()}

            for suffix, kwargs in views.items():
                if dim in (1, 2) and not kwargs.force_view_3d:
                    kwargs.view_2d = True
                    kwargs.scalar_bar_position = [0.04, 0.92, 1.7, 0]
                    if kwargs.grid_vector1 is None:
                        kwargs.grid_vector1 = [1.2, 0, 0]

                    if kwargs.grid_vector2 is None:
                        kwargs.grid_vector2 = [0, -1.2, 0]

                fig_filename = _get_fig_filename(ebase, images_dir, suffix)

                fname = edit_filename(filename, suffix=suffix)
                output('displaying results from "%s"' % fname)
                disp_name = fig_filename.replace(sfepy.data_dir, '')
                output('to "%s"...' % disp_name.lstrip(os.path.sep))

                run_resview_plot(fname, fig_filename, kwargs)

                output('...done')

            remove_files(output_dir)

        output('...done')


def generate_thumbnails(thumbnails_dir, images_dir, scale=0.3):
    """
    Generate thumbnails into `thumbnails_dir` corresponding to images in
    `images_dir`.
    """
    ensure_path(thumbnails_dir + os.path.sep)

    output('generating thumbnails...')
    filenames = glob.glob(os.path.join(images_dir, '*.png'))
    for fig_filename in filenames:
        ebase = fig_filename.replace(sfepy.data_dir, '').lstrip(os.path.sep)
        output('"%s"' % ebase)

        base = os.path.basename(fig_filename)
        thumb_filename = os.path.join(thumbnails_dir, base)

        image.thumbnail(fig_filename, thumb_filename, scale=scale)

    output('...done')


_index = """\
.. _%s-index:

%s
%s

.. toctree::
   :maxdepth: 2

"""

_image = '.. image:: %s'

_include = """\
.. _%s:

%s
%s

**Description**

%s

%s

:download:`source code <%s>`

.. literalinclude:: %s

"""


def generate_rst_files(rst_dir, examples_dir, images_dir, pattern='*.py'):
    """
    Generate Sphinx rst files for examples in `examples_dir` with images
    in `images_dir` and put them into `rst_dir`.

    Returns
    -------
    dir_map : dict
        The directory mapping of examples and corresponding rst files.
    """
    ensure_path(rst_dir + os.path.sep)

    output('generating rst files...')

    dir_map = {}
    for ex_filename in locate_files(pattern, examples_dir):
        if _omit(ex_filename, omits, omit_dirs):
            continue

        ebase = ex_filename.replace(examples_dir, '')[1:]
        base_dir = os.path.dirname(ebase)
        rst_filename = ebase2fbase(ebase) + '.rst'
        dir_map.setdefault(base_dir, []).append((ex_filename, rst_filename))

    for dirname, filenames in dir_map.items():
        filenames = sorted(filenames, key=lambda a: a[1])
        dir_map[dirname] = filenames

    # Main index.
    mfd = open(os.path.join(rst_dir, 'index.rst'), 'w')
    mfd.write(_index % ('examples', 'Examples', '=' * 8))

    for dirname, filenames in ordered_iteritems(dir_map):
        # Subdirectory index.
        ifd = open(os.path.join(rst_dir, '%s-index.rst' % dirname), 'w')
        ifd.write(_index % (dirname + '-examples',
                            dirname, '=' * len(dirname)))

        for ex_filename, rst_filename in filenames:
            full_rst_filename = os.path.join(rst_dir, rst_filename)
            output('"%s"' % rst_filename)
            ebase = ex_filename.replace(examples_dir, '')[1:]

            rst_base = os.path.splitext(rst_filename)[0]

            rst_ex_filename = _make_sphinx_path(ex_filename)
            try:
                docstring = get_default(import_file(ex_filename).__doc__,
                                        'missing description!')

            except KeyboardInterrupt:
                raise

            except:
                output('***** failed! *****')
                docstring = 'example import failed!'

            ifd.write('   %s <%s>\n' % (os.path.basename(ebase), rst_base))
            fig_include = ''
            for fig_filename in _get_fig_filenames(ebase, images_dir):
                rst_fig_filename = _make_sphinx_path(fig_filename)
                if os.path.exists(fig_filename):
                    fig_include += _image % rst_fig_filename + '\n'
                else:
                    output('   warning: figure "%s" not found' % fig_filename)

            # Example rst file.
            fd = open(full_rst_filename, 'w', encoding='utf-8')
            fd.write(_include % (rst_base, ebase, '=' * len(ebase),
                                 docstring,
                                 fig_include,
                                 rst_ex_filename, rst_ex_filename))
            fd.close()

        ifd.close()

        mfd.write('   %s-index\n' % dirname)

    mfd.close()

    output('...done')

    return dir_map


_rst_empty_item = """\
        - ..
"""

_rst_item = """\
      %s - .. figure:: %s
             :target: %s

             :ref:`%s <%s>`
"""

_gallery_table = """\
   .. list-table::
      :align: center
      :class: gallery
"""

_gallery_head = """\
.. only:: html

   .. _gallery-index:

   Gallery
   =======
"""


def generate_gallery(examples_dir, output_filename, doc_dir,
                     rst_dir, thumbnails_dir, dir_map, n_col=3):
    """
    Generate the gallery rst file with thumbnail images and links to
    examples.

    Parameters
    ----------
    output_filename : str
        The output rst file name.
    doc_dir : str
        The top level directory of gallery files.
    rst_dir : str
        The full path to rst files of examples within `doc_dir`.
    thumbnails_dir : str
        The full path to thumbnail images within `doc_dir`.
    dir_map : dict
        The directory mapping returned by `generate_rst_files()`
    n_col : int
        The number of columns in the gallery table.
    """
    output('generating %s...' % output_filename)

    lines = [_gallery_head]

    for dirname, filenames in ordered_iteritems(dir_map):
        title = ['   %s' % dirname.title().replace('_', ' '),
                 '   ' + len(dirname) * '^' + '',
                 _gallery_table]

        llines = []
        icol = 0
        for ex_filename, rst_filename in filenames:
            ebase = ex_filename.replace(examples_dir, '')[1:]
            link = os.path.splitext(rst_filename)[0]

            try:
                thumbnail_filename = next(_get_fig_filenames(ebase,
                                                             thumbnails_dir))
            except StopIteration:
                thumbnail_filename = ''

            if not os.path.isfile(thumbnail_filename):
                # Skip examples with no image (= failed examples).
                output('warning: figure "%s" not found' % thumbnail_filename)
                continue

            thumbnail_name = thumbnail_filename.replace(doc_dir, '..')
            path_to_file = os.path.join(examples_dir, ebase)
            docstring = get_default(import_file(path_to_file).__doc__,
                                    'missing description!')
            docstring = docstring.replace('e.g.', 'eg:')
            docstring = docstring.split('.')
            label = docstring[0].strip()
            label = label.replace('\n', ' ')
            label = label.replace('  ', ' ')

            llines.append(_rst_item % (' ' if icol else '*', thumbnail_name,
                                       link + '.html', label, link))
            icol = (icol + 1) % n_col

        if icol > 0:
            for j in range(icol, n_col):
                llines.append(_rst_empty_item)

        if len(llines) > 0:
            lines += title + llines
        else:
            output('warning: no figures in "%s"' % dirname)

    fd = open(output_filename, 'wt')
    fd.write('\n'.join(lines))
    fd.close()

    output('...done')


helps = {
    'doc_dir': 'top level directory of gallery files',
    'pattern': 'example files search pattern [default: %(default)s]',
    'no_images': 'do not (re)generate images',
    'no_thumbnails': 'do not (re)generate thumbnails',
    'output_filename': 'output file name [default: %(default)s]',
}


def main():
    parser = ArgumentParser(description=__doc__,
                            formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument('--version', action='version', version='%(prog)s')
    parser.add_argument('-d', '--doc-dir', metavar='doc_dir',
                        action='store', dest='doc_dir',
                        default='doc', help=helps['doc_dir'])
    parser.add_argument('-p', '--pattern', metavar='pattern',
                        action='store', dest='pattern',
                        default='*.py', help=helps['pattern'])
    parser.add_argument('-n', '--no-images',
                        action='store_false', dest='images',
                        default=True, help=helps['no_images'])
    parser.add_argument('--no-thumbnails',
                        action='store_false', dest='thumbnails',
                        default=True, help=helps['no_thumbnails'])
    parser.add_argument('-o', '--output', metavar='output_filename',
                        action='store', dest='output_filename',
                        default='gallery.rst',
                        help=helps['output_filename'])
    options = parser.parse_args()

    rst_dir = 'examples'
    doc_dir = os.path.realpath(options.doc_dir)
    examples_dir = os.path.realpath('sfepy/examples')
    full_rst_dir = os.path.join(doc_dir, rst_dir)
    images_dir = os.path.join(doc_dir, 'images/gallery')
    thumbnails_dir = os.path.join(images_dir, 'thumbnails')

    output_filename = os.path.join(full_rst_dir, options.output_filename)

    if options.images:
        generate_images(images_dir, examples_dir, pattern=options.pattern)

    if options.thumbnails:
        generate_thumbnails(thumbnails_dir, images_dir)

    dir_map = generate_rst_files(full_rst_dir, examples_dir, images_dir,
                                 pattern=options.pattern)

    generate_gallery(examples_dir, output_filename, doc_dir,
                     rst_dir, thumbnails_dir, dir_map)


if __name__ == '__main__':
    main()
