r"""
Laplace equation with comments.

Find :math:`t` such that:

.. math::
    \int_{\Omega} c \nabla s \cdot \nabla t
    = 0
    \;, \quad \forall s \;.
"""
#! Poisson Equation
#! ================
#$ \centerline{Example input file, \today}

#! Mesh
#! ----
from sfepy import data_dir

filename_mesh = data_dir + '/meshes/3d/cylinder.mesh'

#! Materials
#! ---------
#$ Here we define just a constant coefficient $c$ of the Poisson equation,
#$ using the 'values' attribute. Other possible attribute is 'function', for
#$ material coefficients computed/obtained at runtime.

material_2 = {
    'name' : 'coef',
    'values' : {'val' : 1.0},
}

#! Regions
#! -------
region_1000 = {
    'name' : 'Omega',
    'select' : 'cells of group 6',
}

region_03 = {
    'name' : 'Gamma_Left',
    'select' : 'vertices in (x < 0.00001)',
    'kind' : 'facet',
}

region_4 = {
    'name' : 'Gamma_Right',
    'select' : 'vertices in (x > 0.099999)',
    'kind' : 'facet',
}

#! Fields
#! ------
#! A field is used mainly to define the approximation on a (sub)domain, i.e. to
#$ define the discrete spaces $V_h$, where we seek the solution.
#!
#! The Poisson equation can be used to compute e.g. a temperature distribution,
#! so let us call our field 'temperature'. On the region 'Omega'
#! it will be approximated using P1 finite elements.

field_1 = {
    'name' : 'temperature',
    'dtype' : 'real',
    'shape' : (1,),
    'region' : 'Omega',
    'approx_order' : 1,
}

#! Variables
#! ---------
#! One field can be used to generate discrete degrees of freedom (DOFs) of
#! several variables. Here the unknown variable (the temperature) is called
#! 't', it's associated DOF name is 't.0' --- this will be referred to
#! in the Dirichlet boundary section (ebc). The corresponding test variable of
#! the weak formulation is called 's'. Notice that the 'dual' item of a test
#! variable must specify the unknown it corresponds to.

variable_1 = {
    'name' : 't',
    'kind' : 'unknown field',
    'field' : 'temperature',
    'order' : 0, # order in the global vector of unknowns
}

variable_2 = {
    'name' : 's',
    'kind' : 'test field',
    'field' : 'temperature',
    'dual' : 't',
}

#! Boundary Conditions
#! -------------------
#! Essential (Dirichlet) boundary conditions can be specified as follows:
ebc_1 = {
    'name' : 't1',
    'region' : 'Gamma_Left',
    'dofs' : {'t.0' : 2.0},
}

ebc_2 = {
    'name' : 't2',
    'region' : 'Gamma_Right',
    'dofs' : {'t.0' : -2.0},
}

#! Equations
#! ---------
#$ The weak formulation of the Poisson equation is:
#$ \begin{center}
#$ Find $t \in V$, such that
#$ $\int_{\Omega} c\ \nabla t : \nabla s = f, \quad \forall s \in V_0$.
#$ \end{center}
#$ The equation below directly corresponds to the discrete version of the
#$ above, namely:
#$ \begin{center}
#$ Find $\bm{t} \in V_h$, such that
#$ $\bm{s}^T (\int_{\Omega_h} c\ \bm{G}^T G) \bm{t} = 0, \quad \forall \bm{s}
#$ \in V_{h0}$,
#$ \end{center}
#$ where $\nabla u \approx \bm{G} \bm{u}$. Below we use $f = 0$ (Laplace
#$ equation).
#! We also define an integral here.
integral_1 = {
    'name' : 'i1',
    'kind' : 'v',
    'order' : 2,
}

equations = {
    'Temperature' : """dw_laplace.i1.Omega( coef.val, s, t ) = 0"""
}

#! Linear solver parameters
#! ---------------------------
#! Use umfpack, if available, otherwise superlu.
solver_0 = {
    'name' : 'ls',
    'kind' : 'ls.scipy_direct',
    'method' : 'auto',
}

#! Nonlinear solver parameters
#! ---------------------------
#! Even linear problems are solved by a nonlinear solver (KISS rule) - only one
#! iteration is needed and the final rezidual is obtained for free.
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

#! Options
#! -------
#! Use them for anything you like... Here we show how to tell which solvers
#! should be used - reference solvers by their names.
options = {
    'nls' : 'newton',
    'ls' : 'ls',
}
