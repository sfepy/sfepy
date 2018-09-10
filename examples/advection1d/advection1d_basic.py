#!/usr/bin/env python
from __future__ import print_function
from __future__ import absolute_import
from argparse import ArgumentParser
import numpy as nm

import sys
sys.path.append('.')

from sfepy.base.base import IndexedStruct
from sfepy.discrete import (FieldVariable, Material, Integral, Function,
                            Equation, Equations, Problem)
from sfepy.discrete.fem import Mesh, FEDomain, Field
from sfepy.terms import Term
from sfepy.discrete.conditions import Conditions, EssentialBC
from sfepy.solvers.ls import ScipyDirect
from sfepy.solvers.nls import Newton
from sfepy.postprocess.viewer import Viewer
from sfepy.mechanics.matcoefs import stiffness_from_lame
from scipy.misc import factorial


mesh1d = Mesh.from_file(r"C:\Users\Lucia Moa\Python_Projects\sfepy\meshes\1d\lspace0-1_10.mesh")
domain1d = FEDomain('domain', mesh1d)
omega1d = domain1d.create_region('Omega', 'all')
gamma11d = domain1d.create_region('Gamma11d',
                                  'vertices in x == 0',
                                  'vertex')
gamma21d = domain1d.create_region('Gamma21d',
                                  'vertices in x > 0.91',
                                  'vertex')
field1d = Field.from_args('fu1d', nm.float64, 'scalar', omega1d,
                          approx_order=2)

u = FieldVariable('u', 'unknown', field1d)
phi = FieldVariable('phi', 'test', field1d, primary_var_name='u')

a = Material('a', val=1)

integral = Integral('i', order=3)

t = Term.new('dw_advect_div_free(a.val, phi, u)', integral, omega1d, a=a, phi=phi, u=u)
eq = Equation('balance', t)
eqs = Equations([eq])
fix_u0 = EssentialBC('fix_u0', gamma11d, {'u.all' : 0.0})
fix_u1 = EssentialBC('fix_u1', gamma11d, {'u.all' : 1.0})
ls = ScipyDirect({})

nls_status = IndexedStruct()
nls = Newton({}, lin_solver=ls, status=nls_status)

pb = Problem('advection1d', equations=eqs)
pb.save_regions_as_groups('regions')

pb.set_bcs(ebcs=Conditions([fix_u0, fix_u1]))
pb.set_solver(nls)
status = IndexedStruct()
state = pb.solve(status=status)
