"""
Configuration classes for acoustic band gaps in a strongly heterogeneous
elastic body.
"""
from __future__ import absolute_import
import numpy as nm

from sfepy.base.base import get_default, import_file, Struct
from sfepy.base.conf import ProblemConf
from sfepy.discrete.fem import MeshIO
import sfepy.discrete.fem.periodic as per
from sfepy.mechanics.matcoefs import stiffness_from_lame, TransformToPlane
from sfepy.homogenization.utils import define_box_regions, get_lattice_volume
import sfepy.homogenization.coefs_base as cb
import sfepy.homogenization.coefs_phononic as cp

per.set_accuracy(1e-8)

def get_pars(dim, lam, mu):
    c = stiffness_from_lame(3, lam, mu)
    if dim == 2:
        tr = TransformToPlane()
        try:
            c = tr.tensor_plane_stress(c3=c)
        except:
            sym = (dim + 1) * dim // 2
            c = nm.zeros((sym, sym), dtype=nm.float64)

    return c

def set_coef_d(variables, ir, ic, mode, pis, corrs_rs):
    mode2var = {'row' : 'u1_m', 'col' : 'u2_m'}

    val = pis.states[ir, ic]['u_m'] + corrs_rs.states[ir, ic]['u_m']

    variables[mode2var[mode]].set_data(val)

class BandGapsConf(Struct):
    """
    Configuration class for acoustic band gaps in a strongly heterogeneous
    elastic body.
    """

    def __init__(self, filename, approx, region_selects, mat_pars, options,
                 evp_options, eigenmomenta_options, band_gaps_options,
                 coefs_save_name='coefs',
                 corrs_save_names=None,
                 incwd=None,
                 output_dir=None, **kwargs):
        Struct.__init__(self, approx=approx, region_selects=region_selects,
                        mat_pars=mat_pars, options=options,
                        evp_options=evp_options,
                        eigenmomenta_options=eigenmomenta_options,
                        band_gaps_options=band_gaps_options,
                        **kwargs)
        self.incwd = get_default(incwd, lambda x: x)

        self.conf = Struct()
        self.conf.filename_mesh = self.incwd(filename)

        output_dir = get_default(output_dir, self.incwd('output'))

        default = {'evp' : 'evp', 'corrs_rs' : 'corrs_rs'}
        self.corrs_save_names = get_default(corrs_save_names,
                                            default)

        io = MeshIO.any_from_filename(self.conf.filename_mesh)
        self.bbox, self.dim = io.read_bounding_box(ret_dim=True)
        rpc_axes = nm.eye(self.dim, dtype=nm.float64) \
                   * (self.bbox[1] - self.bbox[0])

        self.conf.options = options
        self.conf.options.update({
            'output_dir' : output_dir,

            'volume' : {
                'value' : get_lattice_volume(rpc_axes),
            },

            'coefs' : 'coefs',
            'requirements' : 'requirements',

            'coefs_filename' : coefs_save_name,
        })

        self.conf.mat_pars = mat_pars

        self.conf.solvers = self.define_solvers()
        self.conf.regions = self.define_regions()
        self.conf.materials = self.define_materials()
        self.conf.fields = self.define_fields()
        self.conf.variables = self.define_variables()
        (self.conf.ebcs, self.conf.epbcs,
         self.conf.lcbcs, self.all_periodic) = self.define_bcs()
        self.conf.functions = self.define_functions()
        self.conf.integrals = self.define_integrals()

        self.equations, self.expr_coefs = self.define_equations()
        self.conf.coefs = self.define_coefs()
        self.conf.requirements = self.define_requirements()

    def __call__(self):
        return ProblemConf.from_dict(self.conf.__dict__,
                                     import_file(__file__))

    def define_solvers(self):
        solvers = {
            'ls_d' : ('ls.auto_direct', {'use_presolve' : True}),
            'ls_i' : ('ls.scipy_iterative', {
                'method' : 'cg',
                'i_max'      : 1000,
                'eps_a'      : 1e-12,
            }),
            'newton' : ('nls.newton', {
                'i_max' : 1,
                'eps_a' : 1e-4,
            }),
        }

        return solvers

    def define_regions(self):
        regions = {
            'Y' : 'all',
            'Y_m' : self.region_selects.matrix,
            'Y_c' : self.region_selects.inclusion,
            'Gamma_mc': ('r.Y_m *v r.Y_c', 'facet'),
        }

        regions.update(define_box_regions(self.dim,
                                          self.bbox[0], self.bbox[1], 1e-5))

        return regions

    def define_materials(self):
        materials = {
            'm' : ({
                'D_m' : self.mat_pars.D_m,
                'density_m' : self.mat_pars.density_m,
                'D_c' : self.mat_pars.D_c,
                'density_c' : self.mat_pars.density_c,
            }, None, None, {'special_constant' : True}),
        }
        return materials

    def define_fields(self):
        fields = {
            'vector_Y_m' : ('real', self.dim, 'Y_m', self.approx),
            'vector_Y_c' : ('real', self.dim, 'Y_c', self.approx),

            'scalar_Y' : ('real', 1, 'Y', 1),
        }
        return fields

    def define_variables(self):
        variables = {
            'u_m' : ('unknown field', 'vector_Y_m'),
            'v_m' : ('test field', 'vector_Y_m', 'u_m'),
            'Pi'    : ('parameter field', 'vector_Y_m', '(set-to-None)'),
            'u1_m'   : ('parameter field', 'vector_Y_m', '(set-to-None)'),
            'u2_m'   : ('parameter field', 'vector_Y_m', '(set-to-None)'),

            'u_c' : ('unknown field', 'vector_Y_c'),
            'v_c' : ('test field', 'vector_Y_c', 'u_c'),

            'aux'   : ('parameter field', 'scalar_Y', '(set-to-None)'),
        }
        return variables

    def define_bcs(self):
        ebcs = {
            'fixed_corners' : ('Corners', {'u_m.all' : 0.0}),
            'fixed_gamma_mc' : ('Gamma_mc', {'u_c.all' : 0.0}),
        }

        epbcs = {}
        all_periodic = []
        for vn in ['u_m']:
            val = {'%s.all' % vn : '%s.all' % vn}

            epbcs.update({
                'periodic_%s_x' % vn : (['Left', 'Right'], val,
                                        'match_y_line'),
                'periodic_%s_y' % vn : (['Top', 'Bottom'], val,
                                        'match_x_line'),
            })
            all_periodic.extend(['periodic_%s_x' % vn, 'periodic_%s_y' % vn])

        lcbcs = {}

        return ebcs, epbcs, lcbcs, all_periodic

    def define_functions(self):
        functions = {
            'match_x_line' : (per.match_x_line,),
            'match_y_line' : (per.match_y_line,),
        }

        return functions

    def define_integrals(self):
        integrals = {
            'i' : 2,
        }

        return integrals

    def define_equations(self):
        equations = {}
        equations['corrs_rs'] = {
            'balance_of_forces' :
            """dw_lin_elastic.i.Y_m( m.D_m, v_m, u_m )
             = - dw_lin_elastic.i.Y_m( m.D_m, v_m, Pi )""",
        }
        equations['evp'] = {
            'lhs' : """dw_lin_elastic.i.Y_c( m.D_c, v_c, u_c )""",
            'rhs' : """dw_dot.i.Y_c( m.density_c, v_c, u_c )""",
        }

        expr_coefs = {
            'D' : """dw_lin_elastic.i.Y_m( m.D_m, u1_m, u2_m )""",
            'VF' : """ev_volume.i.%s(aux)""",
            'ema' : """ev_integrate.i.Y_c( m.density_c, u_c )""",
        }

        return equations, expr_coefs

    def define_coefs(self):
        from copy import copy

        ema_options = copy(self.eigenmomenta_options)
        ema_options.update({'var_name' : 'u_c'})

        dispersion_options = copy(self.band_gaps_options)
        dispersion_options.update({'log_save_name' : 'dispersion.log'})

        coefs = {
            # Basic.
            'VF' : {
                'regions' : ['Y_m', 'Y_c'],
                'expression' : self.expr_coefs['VF'],
                'class' : cb.VolumeFractions,
            },
            'dv_info' : {
                'requires' : ['c.VF'],
                'region_to_material' : {'Y_m' : ('m', 'density_m'),
                                        'Y_c' : ('m', 'density_c'),},
                'class' : cp.DensityVolumeInfo,
            },

            'eigenmomenta' : {
                'requires' : ['evp', 'c.dv_info'],
                'expression' : self.expr_coefs['ema'],
                'options' : ema_options,
                'class' : cp.Eigenmomenta,
            },
            'M' : {
                'requires' : ['evp', 'c.dv_info', 'c.eigenmomenta'],
                'class' : cp.AcousticMassTensor,
            },
            'band_gaps' : {
                'requires' : ['evp', 'c.eigenmomenta', 'c.M'],
                'options' : self.band_gaps_options,
                'class' : cp.BandGaps,
            },

            # Dispersion.
            'D' : {
                'requires' : ['pis', 'corrs_rs'],
                'expression' : self.expr_coefs['D'],
                'set_variables' : set_coef_d,
                'class' : cb.CoefSymSym,
            },
            'Gamma' : {
                'requires' : ['c.D'],
                'options' : {
                    'mode' : 'simple',
                    'incident_wave_dir' : None,
                },
                'class' : cp.ChristoffelAcousticTensor,
            },
            'dispersion' : {
                'requires' : ['evp', 'c.eigenmomenta', 'c.M', 'c.Gamma'],
                'options' : dispersion_options,
                'class' : cp.BandGaps,
            },
            'polarization_angles' : {
                'requires' : ['c.dispersion'],
                'options' : {
                    'incident_wave_dir' : None,
                },
                'class' : cp.PolarizationAngles,
            },

            # Phase velocity.
            'phase_velocity' : {
                'requires' : ['c.dv_info', 'c.Gamma'],
                'options' : {
                    'eigensolver' : 'eig.sgscipy',
                },
                'class' : cp.PhaseVelocity,
            },
            'filenames' : {},
        }

        return coefs

    def define_requirements(self):
        requirements = {
            # Basic.
            'evp' : {
                'ebcs' : ['fixed_gamma_mc'],
                'epbcs' : None,
                'equations' : self.equations['evp'],
                'save_name' : self.corrs_save_names['evp'],
                'options' : self.evp_options,
                'class' : cp.SimpleEVP,
            },

            # Dispersion.
            'pis' : {
                'variables' : ['u_m'],
                'class' : cb.ShapeDimDim,
            },
            'corrs_rs' : {
                'requires' : ['pis'],
                'ebcs' : ['fixed_corners'],
                'epbcs' : self.all_periodic,
                'equations' : self.equations['corrs_rs'],
                'set_variables' : [('Pi', 'pis', 'u_m')],
                'save_name' : self.corrs_save_names['corrs_rs'],
                'is_linear' : True,
                'class' : cb.CorrDimDim,
                'is_linear' : True,
            },
        }
        return requirements

class BandGapsRigidConf(BandGapsConf):
    """
    Configuration class for acoustic band gaps in a strongly heterogeneous
    elastic body with rigid inclusions.
    """

    def define_regions(self):
        regions = BandGapsConf.define_regions(self)
        regions['Y_cr'] = regions['Y_c']
        regions.update({
            'Y_r' : 'vertices by select_yr',
            'Y_c' : 'r.Y_cr -c r.Y_r',
        })
        return regions

    def define_materials(self):
        materials = BandGapsConf.define_materials(self)
        materials['m'][0].update({
                'D_r' : self.mat_pars.D_r,
                'density_r' : self.mat_pars.density_r,
        })
        return materials

    def define_fields(self):
        fields = {
            'vector_Y_cr' : ('real', self.dim, 'Y_cr', self.approx),

            'scalar_Y' : ('real', 1, 'Y', 1),
        }
        return fields

    def define_variables(self):
        variables = {
            'u' : ('unknown field', 'vector_Y_cr'),
            'v' : ('test field', 'vector_Y_cr', 'u'),

            'aux'   : ('parameter field', 'scalar_Y', '(set-to-None)'),
        }
        return variables

    def define_bcs(self):
        ebcs = {
            'fixed_gamma_mc' : ('Gamma_mc', {'u.all' : 0.0}),
        }
        lcbcs ={
            'rigid' : ('Y_r',{'u.all' : None}, None, 'rigid'),
        }

        return ebcs, {}, lcbcs, []

    def define_functions(self):
        functions = BandGapsConf.define_functions(self)
        functions.update({
            'select_yr' : (self.select_yr,),
        })

        return functions

    def define_equations(self):
        equations = {}

        # dw_lin_elastic.i.Y_r( m.D_r, v, u ) should have no effect!
        equations['evp'] = {
            'lhs' : """dw_lin_elastic.i.Y_c( m.D_c, v, u )
                     + dw_lin_elastic.i.Y_r( m.D_r, v, u )""",
            'rhs' : """dw_dot.i.Y_c( m.density_c, v, u )
                     + dw_dot.i.Y_r( m.density_r, v, u )""",
        }

        expr_coefs = {
            'VF' : """ev_volume.i.%s(aux)""",
            'ema' : """ev_integrate.i.Y_c( m.density_c, u )
                     + ev_integrate.i.Y_r( m.density_r, u )""",
        }

        return equations, expr_coefs

    def define_coefs(self):
        from copy import copy

        ema_options = copy(self.eigenmomenta_options)
        ema_options.update({'var_name' : 'u'})

        coefs = {
            # Basic.
            'VF' : {
                'regions' : ['Y_m', 'Y_cr', 'Y_c', 'Y_r'],
                'expression' : self.expr_coefs['VF'],
                'class' : cb.VolumeFractions,
            },
            'dv_info' : {
                'requires' : ['c.VF'],
                'region_to_material' : {'Y_m' : ('m', 'density_m'),
                                        'Y_c' : ('m', 'density_c'),
                                        'Y_r' : ('m', 'density_r'),},
                'class' : cp.DensityVolumeInfo,
            },

            'eigenmomenta' : {
                'requires' : ['evp', 'c.dv_info'],
                'expression' : self.expr_coefs['ema'],
                'options' : ema_options,
                'class' : cp.Eigenmomenta,
            },
            'M' : {
                'requires' : ['evp', 'c.dv_info', 'c.eigenmomenta'],
                'class' : cp.AcousticMassTensor,
            },
            'band_gaps' : {
                'requires' : ['evp', 'c.eigenmomenta', 'c.M'],
                'options' : self.band_gaps_options,
                'class' : cp.BandGaps,
            },

            'filenames' : {},
        }

        return coefs

    def define_requirements(self):
        requirements = {
            # Basic.
            'evp' : {
                'ebcs' : ['fixed_gamma_mc'],
                'epbcs' : None,
                'lcbcs' : ['rigid'],
                'equations' : self.equations['evp'],
                'save_name' : self.corrs_save_names['evp'],
                'options' : self.evp_options,
                'class' : cp.SimpleEVP,
            },
        }
        return requirements


def clip(data, plot_range):
    return nm.clip(data, *plot_range)

def clip_sqrt(data, plot_range):
    return nm.clip(nm.sqrt(data), *plot_range)

def normalize(data, plot_range):
    aux = nm.arctan(data)
    return clip(aux, plot_range)
