# -*- coding: utf-8 -*-
# last revision: 25.02.2008
filename_mesh = 'database/simple.vtk'

###############################
#
# This file models a cylinder that is fixed at both ends. The input to the
# problem is a specified displacement of 0.01 in the z direction for points in
# the region labeled SomewhereTop.  This boundary condition is named
# PerturbedSurface.  The region SomewhereTop is specified as those nodes for
# which
#              (z > 0.017) & (x > 0.01) & (x <  0.08).
# The output is the displacement for each node, saved by default to
# simple_out.vtk. The material is linear elastic and its properties are
# specified as Lamé parameters (see
# http://en.wikipedia.org/wiki/Lam%C3%A9_parameters)
#
###############################

field_1 = {
    'name' : '3_displacement',
    'dim' : (3,1),
    'domain' : 'Y',
    'bases' : {'Y' : '3_4_P1'}
}

material_1 = {
    'name' : 'solid',
    'mode' : 'here',
    'region' : 'Y',
    'lame' : {'lambda' : 1e1, 'mu' : 1e0},
}

variable_1 = {
    'name' : 'u',
    'kind' : 'unknown field',
    'field' : '3_displacement',
    'order' : 0,
}
variable_2 = {
    'name' : 'v',
    'kind' : 'test field',
    'field' : '3_displacement',
    'dual' : 'u',
}

# Whole domain $Y$.
region_1000 = {
    'name' : 'Y',
    'select' : 'all',
}

# EBC regions.
region_1 = {
    'name' : 'Left',
    'select' : 'nodes in (x < 0.001)'
}
region_2 = {
    'name' : 'Right',
    'select' : 'nodes in (x > 0.099)'
}
region_3 = {
    'name' : 'SomewhereTop',
    'select' : 'nodes in (z > 0.017) & (x > 0.01) & (x < 0.08)'
}

ebc_1 = {
    'name' : 'Left',
    'region' : 'Left',
    'dofs' : {'u.all' : 0.0},
}
ebc_2 = {
    'name' : 'Right',
    'region' : 'Right',
    'dofs' : {'u.all' : 0.0},
}
ebc_3 = {
    'name' : 'PerturbedSurface',
    'region' : 'SomewhereTop',
    'dofs' : {'u.2' : 0.01},
}

integral_1 = {
    'name' : 'i1',
    'kind' : 'v',
    'quadrature' : 'gauss_o1_d3',
}
equations = {
    'balance_of_forces' : """dw_lin_elastic_iso.i1.Y( solid.lame, v, u ) = 0""",
}

solver_0 = {
    'name' : 'ls',
    'kind' : 'ls.umfpack',
}

solver_1 = {
    'name' : 'newton',
    'kind' : 'nls.newton',

    'i_max'      : 1,
    'eps_a'      : 1e-10,
    'eps_r'      : 1.0,
    'macheps'   : 1e-16,
    'lin_red'    : 1e-2, # Linear system error < (eps_a * lin_red).
    'ls_red'     : 0.1,
    'ls_red_warp' : 0.001,
    'ls_on'      : 1.1,
    'ls_min'     : 1e-5,
    'check'     : 0,
    'delta'     : 1e-6,
    'is_plot'    : False,
    'problem'   : 'nonlinear', # 'nonlinear' or 'linear' (ignore i_max)
}

##
# FE assembling parameters.
fe = {
    'chunk_size' : 1000
}

##
# Functions.
from valec import *

##
# Make 'cinc' refer to a cinc_* function according to the mesh file name.
import os.path as op
trunk = op.splitext( op.basename( filename_mesh ) )[0]
cinc = eval( 'cinc_' + trunk )
del op, trunk
