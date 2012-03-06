r"""
Stabilized Navier-Stokes problem with grad-div, SUPG and PSPG stabilization
solved by a custom Oseen solver, see [1].

[1] G. Matthies and G. Lube. On streamline-diffusion methods of inf-sup stable
discretisations of the generalised Oseen problem. Number 2007-02 in Preprint
Series of Institut fuer Numerische und Angewandte Mathematik,
Georg-August-Universitaet Goettingen, 2007.

Find :math:`\ul{u}`, :math:`p` such that:

.. math::
    \begin{array}{l}
    \int_{\Omega} \nu\ \nabla \ul{v} : \nabla \ul{u}
    \int_{\Omega} ((\ul{b} \cdot \nabla) \ul{u}) \cdot \ul{v}
    - \int_{\Omega} p\ \nabla \cdot \ul{v} \\
    + \gamma \int_{\Omega} (\nabla\cdot\ul{u}) \cdot (\nabla\cdot\ul{v}) \\
    + \sum_{K \in \Ical_h}\int_{T_K} \delta_K\ ((\ul{b} \cdot \nabla)
      \ul{u})\cdot ((\ul{b} \cdot \nabla) \ul{v}) \\
    + \sum_{K \in \Ical_h}\int_{T_K} \delta_K\ \nabla p\cdot ((\ul{b} \cdot
      \nabla) \ul{v})
    = 0
    \;, \quad \forall \ul{v} \;,
    \end{array}

    \begin{array}{l}
    \int_{\Omega} q\ \nabla \cdot \ul{u} \\
    + \sum_{K \in \Ical_h}\int_{T_K} \tau_K\ ((\ul{b} \cdot \nabla) \ul{u})
      \cdot \nabla q \\
    + \sum_{K \in \Ical_h}\int_{T_K} \tau_K\ \nabla p \cdot \nabla q
    = 0
    \;, \quad \forall q \;.
    \end{array}
"""
from sfepy.base.base import get_default, Struct
from sfepy import data_dir

filename_mesh = data_dir + '/meshes/3d/elbow2.mesh'

options = {
    'solution' : 'steady',
    'nls' : 'oseen',
    'ls' : 'ls',
}

regions = {
    'Omega' : ('all', {}),
    'Walls' : ('nodes of surface -n (r.Outlet +n r.Inlet)',
               {'can_cells' : False}),
    'Inlet' : ('nodes by cinc0', {'can_cells' : False}),
    'Outlet' : ('nodes by cinc1', {'can_cells' : False}),
}

fields = {
    'velocity' : ('real', 3, 'Omega', 1),
    'pressure' : ('real', 1, 'Omega', 1),
}

variables = {
    'u'   : ('unknown field',   'velocity', 0),
    'v'   : ('test field',      'velocity', 'u'),
    'b'   : ('parameter field', 'velocity', 'u'),
    'p'   : ('unknown field',   'pressure', 1),
    'q'   : ('test field',      'pressure', 'p'),
}

ebcs = {
    'Walls_velocity' : ('Walls', {'u.all' : 0.0}),
    'Inlet_velocity' : ('Inlet', {'u.1' : 1.0, 'u.[0,2]' : 0.0}),
}

materials = {
    'fluid' : ({'viscosity' : 1.25e-5,
                'density' : 1e0},),
    'stabil' : 'stabil',
}

integrals = {
    'i1' : ('v', 'gauss_o2_d3'),
    'i2' : ('v', 'gauss_o3_d3'),
}

##
# Stationary Navier-Stokes equations with grad-div, SUPG and PSPG stabilization.
equations = {
    'balance' :
    """  dw_div_grad.i2.Omega( fluid.viscosity, v, u )
       + dw_lin_convect.i2.Omega( v, b, u )
       - dw_stokes.i1.Omega( v, p )
       + dw_st_grad_div.i1.Omega( stabil.gamma, v, u )
       + dw_st_supg_c.i1.Omega( stabil.delta, v, b, u )
       + dw_st_supg_p.i1.Omega( stabil.delta, v, b, p )
       = 0""",
    'incompressibility' :
    """  dw_stokes.i1.Omega( u, q )
       + dw_st_pspg_c.i1.Omega( stabil.tau, q, b, u )
       + dw_st_pspg_p.i1.Omega( stabil.tau, q, p )
       = 0""",
}

solver_1 = {
    'name' : 'oseen',
    'kind' : 'nls.oseen',

    'needs_problem_instance' : True,
    'stabil_mat' : 'stabil',

    'adimensionalize' : False,
    'check_navier_stokes_rezidual' : False,

    'i_max'      : 10,
    'eps_a'      : 1e-8,
    'eps_r'      : 1.0,
    'macheps'    : 1e-16,
    'lin_red'    : 1e-2, # Linear system error < (eps_a * lin_red).
    'is_plot'    : False,

    # Uncomment the following to get a convergence log.
    ## 'log'        : {'text' : 'oseen_log.txt',
    ##                 'plot' : 'oseen_log.png'},
}

solver_2 = {
    'name' : 'ls',
    'kind' : 'ls.scipy_direct',
}

##
# Functions.
import os.path as op
import sys

import numpy as nm

sys.path.append(data_dir) # Make installed example work.
import examples.navier_stokes.utils as utils

cinc_name = 'cinc_' + op.splitext(op.basename(filename_mesh))[0]
cinc = getattr(utils, cinc_name)

class StabilizationFunction(Struct):

    gamma = None
    delta = None
    tau = None
    tau_red = 1.0e-0 # <= 1.0; if tau is None: tau = tau_red * delta
    tau_mul = 1.0
    delta_mul = 1.0e-0
    gamma_mul = 1.0e0
    # 'edge': longest edge 'volume': volume-based, 'max': max. of
    # previous
    diameter_mode = 'max' # 'edge', 'volume', 'max'

    def setup(self, problem):
        """
        Setup common problem-dependent data.
        """
        variables = problem.get_variables()

        # Map names for Oseen solver.
        self.ns = {'p' : 'p', 'q' : 'q',
                   'u' : 'u', 'b' : 'b', 'v' : 'v',
                   'fluid' : 'fluid', 'omega' : 'omega',
                   'i1' : 'i1', 'i2' : 'i2'}

        # Indices to the state vector.
        ii = {}
        ii['u'] = variables.get_indx('u')
        ii['us'] = variables.get_indx('u', stripped=True)
        ii['ps'] = variables.get_indx('p', stripped=True)
        self.ii = ii

        materials = problem.get_materials()

        # The viscosity.
        fluid_mat = materials['fluid']
        self.viscosity = fluid_mat.function()['viscosity']

        # The Friedrich's constant.
        self.c_friedrichs = problem.domain.get_diameter()
        self.sigma = 1e-12 # 1 / dt.

        # Element diameter modes.
        dm = {'edge': 0, 'volume': 1, 'max': 2}[self.diameter_mode]

        self.b_norm = 1.0

        term = problem.equations['balance'].terms['dw_lin_convect']
        region = term.region
        var = variables['u']
        diameters2 = []
        for ig in term.iter_groups():
            vg, _ = term.get_mapping(var)
            cells = region.get_cells(ig)
            d2 = problem.domain.get_element_diameters(ig, cells, vg, dm)
            diameters2.append(d2)

        self.diameters2 = nm.concatenate(diameters2)

    def get_maps(self):
        return self.ns, self.ii

    def __call__(self, ts, coor, mode=None, term=None, problem=None,
                 b_norm=None, **kwargs):
        if mode != 'qp': return

        if not hasattr(self, 'viscosity'):
            self.setup(problem)

        # Update stored b_norm.
        self.b_norm = get_default(b_norm, self.b_norm)

        print '|b|_max (mat_fun):', self.b_norm
        gamma = self.viscosity + self.b_norm * self.c_friedrichs

        data = {}
        if self.gamma is None:
            data['gamma'] = self.gamma_mul * gamma

        else:
            data['gamma'] = nm.asarray(self.gamma_mul * self.gamma,
                                       dtype=nm.float64)
        data['gamma'] = nm.tile(data['gamma'], (coor.shape[0], 1, 1))

        if self.delta is None:
            val1 = min(1.0, 1.0 / self.sigma)
            val2 = self.sigma * self.c_friedrichs**2
            val3 = (self.b_norm**2) \
                   * min((self.c_friedrichs**2) / self.viscosity,
                         1.0 / self.sigma)

            n_qp = coor.shape[0] / self.diameters2.shape[0]
            diameters2 = nm.repeat(self.diameters2, n_qp)
            diameters2.shape = diameters2.shape + (1, 1)

            data['delta'] = self.delta_mul * val1 \
                            * diameters2 / (data['gamma'] + val2 + val3)

        else:
            val = nm.asarray(self.delta_mul * self.delta, dtype=nm.float64)
            data['delta'] = nm.tile(val, (coor.shape[0], 1, 1))

        if self.tau is None:
            data['tau'] = self.tau_red * data['delta']

        else:
            data['tau'] = nm.asarray(self.tau_mul * self.tau, dtype=nm.float64)
            data['tau'] = nm.tile(data['tau'], (coor.shape[0], 1, 1))

        return data

functions = {
    'cinc0' : (lambda coors, domain=None: cinc(coors, 0),),
    'cinc1' : (lambda coors, domain=None: cinc(coors, 1),),
    'stabil' : (StabilizationFunction(),),
}
