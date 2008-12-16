from sfepy.fem.periodic import *
from sfepy.homogenization.coefs import CorrDimDim, ShapeDimDim, CoefSymSym

from regions import get_box_regions

def expand_regions( eqs, expand ):
    out = {}
    for key, val in eqs.iteritems():
        out[key] = val % expand
    return out

class CorrectorsRS( CorrDimDim ):
    """Steady state correctors $\bar{\omega}^{rs}$."""

    def get_variables( self, ir, ic, data ):
        """data: pis"""
        pis = data[self.requires[0]]
        yield (self.variables[-1], pis[ir,ic])

class ElasticCoef( CoefSymSym ):
    """Homogenized elastic tensor $E_{ijkl}$."""

    mode2var = {'row' : 0, 'col' : 1}

    def get_variables( self, problem, ir, ic, data, mode ):
        pis, corrs = [data[ii] for ii in self.requires]

        var_name = self.variables[self.mode2var[mode]]
        u_name = problem.variables[var_name].primary_var_name

        indx = corrs.di.indx[u_name]
        omega = corrs.states[ir,ic][indx]
        pi = pis[ir,ic] + omega
        yield (var_name, pi)

expr_elastic = """dw_lin_elastic.i2.%s( matrix.D, Pi1, Pi2 )"""

eq_rs = {
    'eq' : """dw_lin_elastic.i2.%s( matrix.D, v1, u1 )
            + dw_lin_elastic.i2.%s( matrix.D, v1, Pi ) = 0""",
}

def define_input( filename, region, dim, geom ):
    """Uses materials, fe of master file, merges regions."""
    filename_mesh = filename
    
    options = {
        'coefs' : 'coefs',
        'requirements' : 'requirements',
    }

    coefs = {
        'elastic' : {
            'requires' : ['pis', 'corrs_phono_rs'],
            'variables' : ['Pi1', 'Pi2', 'pi1', 'pi2'],
            'expression' : expr_elastic % region,
            'class' : ElasticCoef,
        },
    }

    all_periodic = ['periodic_%s' % ii for ii in ['x', 'y', 'z'][:dim] ]

    requirements = {
        'pis' : {
            'variables' : ['u1'],
            'class' : ShapeDimDim,
        },
        'corrs_phono_rs' : {
            'requires' : ['pis'],
            'variables' : ['u1', 'v1', 'Pi'],
            'ebcs' : ['fixed_u'],
            'epbcs' : all_periodic,
            'equations' : expand_regions( eq_rs, (region, region) ),
            'class' : CorrectorsRS,
            'save_name' : 'corrs_phono',
        },
    }

    integral_1 = {
        'name' : 'i2',
        'kind' : 'v',
        'quadrature' : 'gauss_o3_d%d' % dim,
    }

    field_10 = {
        'name' : 'displacement_matrix',
        'dim' : (dim,1),
        'domain' : region,
        'bases' : {region : '%s_P1' % geom}
    }

    variables = {
        'u1' : ('unknown field', 'displacement_matrix', 0),
        'v1' : ('test field', 'displacement_matrix', 'u1'),
        'Pi' : ('parameter field', 'displacement_matrix', 'u1'),
        'Pi1' : ('parameter field', 'displacement_matrix', 'u1'),
        'Pi2' : ('parameter field', 'displacement_matrix', 'u1'),
    }

    if filename.find( 'mesh_circ21' ) >= 0:
        sizes = (0.499, 0.499)
    elif filename.find( 'cube_cylinder_centered' ) >= 0:
        sizes = (0.499, 0.499, 0.499)

    regions = get_box_regions( dim, sizes )

    ebcs = {
        'fixed_u' : ('Corners', {'u1.all' : 0.0}),
    }

    ##
    # Periodic boundary conditions.
    if dim == 3:
        epbc_10 = {
            'name' : 'periodic_x',
            'region' : ['Left', 'Right'],
            'dofs' : {'u1.all' : 'u1.all'},
            'match' : 'match_x_plane',
        }
        epbc_11 = {
            'name' : 'periodic_y',
            'region' : ['Near', 'Far'],
            'dofs' : {'u1.all' : 'u1.all'},
            'match' : 'match_y_plane',
        }
        epbc_12 = {
            'name' : 'periodic_z',
            'region' : ['Top', 'Bottom'],
            'dofs' : {'u1.all' : 'u1.all'},
            'match' : 'match_z_plane',
        }
    else:
        epbc_10 = {
            'name' : 'periodic_x',
            'region' : ['Left', 'Right'],
            'dofs' : {'u1.all' : 'u1.all'},
            'match' : 'match_y_line',
        }
        epbc_11 = {
            'name' : 'periodic_y',
            'region' : ['Top', 'Bottom'],
            'dofs' : {'u1.all' : 'u1.all'},
            'match' : 'match_x_line',
        }

    solver_0 = {
        'name' : 'ls',
        'kind' : 'ls.umfpack', # Direct solver.
    }

    solver_1 = {
        'name' : 'newton',
        'kind' : 'nls.newton',

        'i_max'      : 2,
    }

    return locals()
