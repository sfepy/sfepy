#!/usr/bin/env python
"""
Compute homogenized elastic coefficients for a given microstructure.
"""
from __future__ import print_function
from __future__ import absolute_import
from argparse import ArgumentParser
import sys
import six
sys.path.append('.')

import numpy as nm

from sfepy import data_dir
import sfepy.discrete.fem.periodic as per
from sfepy.homogenization.utils import define_box_regions

def define_regions(filename):
    """
    Define various subdomains for a given mesh file.
    """
    regions = {}
    dim = 2

    regions['Y'] = 'all'

    eog = 'cells of group %d'
    if filename.find('osteonT1') >= 0:
        mat_ids = [11, 39, 6, 8, 27, 28, 9, 2, 4, 14, 12, 17, 45, 28, 15]
        regions['Ym'] = ' +c '.join((eog % im) for im in  mat_ids)
        wx = 0.865
        wy = 0.499

    regions['Yc'] = 'r.Y -c r.Ym'

    # Sides and corners.
    regions.update(define_box_regions(2, (wx, wy)))

    return dim, regions

def get_pars(ts, coor, mode=None, term=None, **kwargs):
    """
    Define material parameters: :math:`D_ijkl` (elasticity), in a given region.
    """
    if mode == 'qp':
        dim = coor.shape[1]
        sym = (dim + 1) * dim // 2

        out = {}

        # in 1e+10 [Pa]
        lam = 1.7
        mu = 0.3
        o = nm.array([1.] * dim + [0.] * (sym - dim), dtype = nm.float64)
        oot = nm.outer(o, o)
        out['D'] = lam * oot + mu * nm.diag(o + 1.0)

        for key, val in six.iteritems(out):
            out[key] = nm.tile(val, (coor.shape[0], 1, 1))

        channels_cells = term.region.domain.regions['Yc'].cells
        n_cell = term.region.get_n_cells()
        val = out['D'].reshape((n_cell, -1, 3, 3))
        val[channels_cells] *= 1e-1

        return out

##
# Mesh file.
filename_mesh = data_dir + '/meshes/2d/special/osteonT1_11.mesh'

##
# Define regions (subdomains, boundaries) - $Y$, $Y_i$, ...
# depending on a mesh used.
dim, regions = define_regions(filename_mesh)

functions = {
    'get_pars' : (lambda ts, coors, **kwargs:
                  get_pars(ts, coors, **kwargs),),
    'match_x_plane' : (per.match_x_plane,),
    'match_y_plane' : (per.match_y_plane,),
    'match_z_plane' : (per.match_z_plane,),
    'match_x_line' : (per.match_x_line,),
    'match_y_line' : (per.match_y_line,),
}

##
# Define fields: 'displacement' in $Y$,
# 'pressure_m' in $Y_m$.
fields = {
    'displacement' : ('real', dim, 'Y', 1),
}

##
# Define corrector variables: unknown displaements: uc, test: vc
# displacement-like variables: Pi, Pi1, Pi2
variables = {
    'uc'       : ('unknown field',   'displacement', 0),
    'vc'       : ('test field',      'displacement', 'uc'),
    'Pi'       : ('parameter field', 'displacement', 'uc'),
    'Pi1'      : ('parameter field', 'displacement', None),
    'Pi2'      : ('parameter field', 'displacement', None),
}

##
# Periodic boundary conditions.
if dim == 3:
    epbcs = {
        'periodic_x' : (['Left', 'Right'], {'uc.all' : 'uc.all'},
                        'match_x_plane'),
        'periodic_y' : (['Near', 'Far'], {'uc.all' : 'uc.all'},
                        'match_y_plane'),
        'periodic_z' : (['Top', 'Bottom'], {'uc.all' : 'uc.all'},
                        'match_z_plane'),
    }
else:
    epbcs = {
        'periodic_x' : (['Left', 'Right'], {'uc.all' : 'uc.all'},
                        'match_y_line'),
        'periodic_y' : (['Bottom', 'Top'], {'uc.all' : 'uc.all'},
                        'match_x_line'),
    }

##
# Dirichlet boundary conditions.
ebcs = {
    'fixed_u' : ('Corners', {'uc.all' : 0.0}),
}

##
# Material defining constitutive parameters of the microproblem.
materials = {
    'm' : 'get_pars',
}

##
# Numerical quadratures for volume (i3 - order 3) integral terms.
integrals = {
    'i3' : 3,
}

##
# Homogenized coefficients to compute.
def set_elastic(variables, ir, ic, mode, pis, corrs_rs):
    mode2var = {'row' : 'Pi1', 'col' : 'Pi2'}

    val = pis.states[ir, ic]['uc'] + corrs_rs.states[ir, ic]['uc']

    variables[mode2var[mode]].set_data(val)

coefs = {
    'E' : {
        'requires' : ['pis', 'corrs_rs'],
        'expression' : 'dw_lin_elastic.i3.Y(m.D, Pi1, Pi2)',
        'set_variables' : set_elastic,
    },
}

all_periodic = ['periodic_%s' % ii for ii in ['x', 'y', 'z'][:dim] ]
requirements = {
    'pis' : {
        'variables' : ['uc'],
    },
    ##
    # Steady state correctors $\bar{\omega}^{rs}$.
    'corrs_rs' : {
        'requires' : ['pis'],
        'save_variables' : ['uc'],
        'ebcs' : ['fixed_u'],
        'epbcs' : all_periodic,
        'equations' : {'eq' : """dw_lin_elastic.i3.Y(m.D, vc, uc)
                             = - dw_lin_elastic.i3.Y(m.D, vc, Pi)"""},
        'set_variables' : [('Pi', 'pis', 'uc')],
        'save_name' : 'corrs_elastic',
        'is_linear' : True,
    },
}

##
# Solvers.
solvers = {
    'ls' : ('ls.auto_direct', {'use_presolve' : True}),
    'newton' : ('nls.newton', {
        'i_max' : 1,
        'eps_a' : 1e-8,
        'eps_r' : 1e-2,
    })
}

############################################
# Mini-application below, computing the homogenized elastic coefficients.
helps = {
    'no_pauses' : 'do not make pauses',
}

def main():
    import os
    from sfepy.base.base import spause, output
    from sfepy.base.conf import ProblemConf, get_standard_keywords
    from sfepy.discrete import Problem
    import sfepy.homogenization.coefs_base as cb

    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--version', action='version', version='%(prog)s')
    parser.add_argument('-n', '--no-pauses',
                        action="store_true", dest='no_pauses',
                        default=False, help=helps['no_pauses'])
    options = parser.parse_args()

    if options.no_pauses:
        def spause(*args):
            output(*args)

    nm.set_printoptions(precision=3)

    spause(r""">>>
First, this file will be read in place of an input
(problem description) file.
Press 'q' to quit the example, press any other key to continue...""")
    required, other = get_standard_keywords()
    required.remove('equations')
    # Use this file as the input file.
    conf = ProblemConf.from_file(__file__, required, other)
    print(list(conf.to_dict().keys()))
    spause(r""">>>
...the read input as a dict (keys only for brevity).
['q'/other key to quit/continue...]""")

    spause(r""">>>
Now the input will be used to create a Problem instance.
['q'/other key to quit/continue...]""")
    problem = Problem.from_conf(conf, init_equations=False)
    # The homogenization mini-apps need the output_dir.
    output_dir = ''
    problem.output_dir = output_dir
    print(problem)
    spause(r""">>>
...the Problem instance.
['q'/other key to quit/continue...]""")

    spause(r""">>>
The homogenized elastic coefficient $E_{ijkl}$ is expressed
using $\Pi$ operators, computed now. In fact, those operators are permuted
coordinates of the mesh nodes.
['q'/other key to quit/continue...]""")
    req = conf.requirements['pis']
    mini_app = cb.ShapeDimDim('pis', problem, req)
    mini_app.setup_output(save_formats=['vtk'])
    pis = mini_app()
    print(pis)
    spause(r""">>>
...the $\Pi$ operators.
['q'/other key to quit/continue...]""")

    spause(r""">>>
Next, $E_{ijkl}$ needs so called steady state correctors $\bar{\omega}^{rs}$,
computed now.
['q'/other key to quit/continue...]""")
    req = conf.requirements['corrs_rs']

    save_name = req.get('save_name', '')
    name = os.path.join(output_dir, save_name)

    mini_app = cb.CorrDimDim('steady rs correctors', problem, req)
    mini_app.setup_output(save_formats=['vtk'])
    corrs_rs = mini_app(data={'pis': pis})
    print(corrs_rs)
    spause(r""">>>
...the $\bar{\omega}^{rs}$ correctors.
The results are saved in: %s.%s

Try to display them with:

   sfepy-view %s.%s

['q'/other key to quit/continue...]""" % (2 * (name, problem.output_format)))

    spause(r""">>>
Then the volume of the domain is needed.
['q'/other key to quit/continue...]""")
    volume = problem.evaluate('ev_volume.i3.Y(uc)')
    print(volume)

    spause(r""">>>
...the volume.
['q'/other key to quit/continue...]""")

    spause(r""">>>
Finally, $E_{ijkl}$ can be computed.
['q'/other key to quit/continue...]""")
    mini_app = cb.CoefSymSym('homogenized elastic tensor',
                             problem, conf.coefs['E'])
    c_e = mini_app(volume, data={'pis': pis, 'corrs_rs' : corrs_rs})
    print(r""">>>
The homogenized elastic coefficient $E_{ijkl}$, symmetric storage
with rows, columns in 11, 22, 12 ordering:""")
    print(c_e)

if __name__ == '__main__':
    main()
