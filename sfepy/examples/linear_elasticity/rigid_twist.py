r"""
Linear elasticity example demonstrating use of the 'rigid_twist' LCBCs.

Find :math:`\ul{u}` such that:

.. math::
    \int_{\Omega} D_{ijkl}\ e_{ij}(\ul{v}) e_{kl}(\ul{u})
    = 0
    \;, \quad \forall \ul{v} \;,

and displacements in ``'RigidM'`` and ``'RigidS'`` regions are tied using the
:class:`RigidTwistOperator
<sfepy.discrete.fem.lcbc_operators.RigidTwistOperator>`.

Usage Examples
--------------

- Run with the default parameters::

    sfepy-run sfepy/examples/linear_elasticity/rigid_twist.py
    sfepy-view output/rigid_twist/cut-cylinder.*.vtk -f u:wu:f1:p0 1:vw:wu:f1:p0 1:vw:p0

- Use compression instead of tension, decrease top plane shift, tweak thread
  height parameter::

    sfepy-run sfepy/examples/linear_elasticity/rigid_twist.py -d solver=auto,shift=-0.01,thread=0.1,refine=0
    sfepy-view output/rigid_twist/cut-cylinder.*.vtk -f u:wu:f10:p0 1:vw:wu:f10:p0 1:vw:p0
"""
import numpy as nm

from sfepy import data_dir
from sfepy.mechanics.matcoefs import stiffness_from_youngpoisson

def define(
        E=5e6,
        nu=0.45,
        thread=0.34,
        shift=0.05,
        refine=0,
        solver='auto',
        output_dir='output/rigid_twist',
        **kwargs,
):
    filename_mesh = data_dir + '/meshes/3d/cut-cylinder.vtk'

    cm, cs = -0.025, 0.025
    eps = 1e-8
    cm0, cm1 = cm - eps, cm + eps
    cs0, cs1 = cs - eps, cs + eps

    centre = [0, 0, 0]
    c = 'z'
    axis = [0.0, 0.0, 1.0]

    options = {
        'nls' : 'newton',
        'ls' : solver,
        'refinement_level' : refine,
        'output_dir': output_dir,
    }

    regions = {
        'Omega' :  'all',
        'Bottom' : (f'vertices in ({c} < -0.499999)', 'facet'),
        'Top' : (f'vertices in ({c} > 0.499999)', 'facet'),
        'RigidM' : (f'vertices in ({c} > {cm0}) & ({c} < {cm1})', 'facet'),
        'RigidS' : (f'vertices in ({c} > {cs0}) & ({c} < {cs1})', 'facet'),
    }

    fields = {
        'fu' : ('real', 'vector', 'Omega', 1),
    }

    variables = {
        'u' : ('unknown field', 'fu', 3),
        'v' : ('test field', 'fu', 'u'),
    }

    materials = {
        'm' : ({
            'D' : stiffness_from_youngpoisson(3, young=E, poisson=nu),
        },)
    }

    ebcs = {
        'bottom': ('Bottom', {'u.all': 0.0}),
        'top': ('Top', {'u.all': 'move_top'}),
    }

    def move_top(ts, coors, bc, problem, **kwargs):
        val = nm.zeros_like(coors)
        ic = dict(x=0, y=1, z=2)[c]
        val[:, ic] = ts.nt * shift
        return val

    functions = {
        'move_top' : (move_top,),
    }

    lcbcs = {
        'rt' : ('RigidM', {'u.all' : None}, None, 'rigid_twist',
                'RigidS', centre, axis, thread),
    }

    integrals = {
        'i': 2,
    }

    equations = {
        'balance_of_forces' :
        """dw_lin_elastic.i.Omega(m.D, v, u) = 0""",
    }

    solvers = {
        'pypardiso': ('ls.pypardiso', {}),
        'auto': ('ls.auto_direct', {
            'use_presolve' : True,
        }),
        'newton' : ('nls.newton', {
            'i_max' : 1,
            'eps_a' : 1e-8,
            'lin_red' : None,
            'report_status' : True,
            'is_linear' : True,
        }),
        'ts' : ('ts.simple', {
            't0'     : 0.0,
            't1'     : 1.0,
            'dt'     : None,
            'n_step' : 5,
            'quasistatic' : True,
            'verbose' : 1,
        }),
    }

    return locals()
