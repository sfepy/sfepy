"""
Solves Poisson equation. Demonstrates parametric study capabilities of
Application classes.

In particular (in pseudo-LaTeX):

c \nabla t = f in \Omega,
t = 2 on \Gamma_1,  t = -2 on \Gamma_2,
f = 1 in \Omega_1, f = 0 otherwise,

\Omega is a square domain, \Omega_1 \in \Omega is a cicular domain.
Now let's see what happens if \Omega_1 diameter changes.

Run
$ ./simple.py <this file>
and then look in 'output/r_omega1' directory.

Remark: this simple case could be achieved also by defining \Omega_1 by a
time-dependent function and solve the static problem as a time-dependent
problem. However, the approach below is much more general.
"""
import os
import numpy as nm
from sfepy.base.base import output, set_output_prefix, pause, debug

# Mesh.
filename_mesh = 'database/square_circ.vtk'

# Options. The value of 'parametric_hook' is the function that does the
# parametric study.
options = {
    'nls' : 'newton', # Nonlinear solver
    'ls' : 'ls', # Linear solver

    'parametric_hook' : 'vary_omega1_size',
    'output_dir' : 'output/r_omega1',
}

# Domain and subdomains.
default_diameter = 0.25
regions = {
    'Omega' : ('all', {}),
    'Gamma_1' : ('nodes in (x < -0.999)', {}),
    'Gamma_2' : ('nodes in (x > 0.999)', {}),
    'Omega_1' : ('nodes by select_circ( x, y, z, %f )' % default_diameter, {}),
}

# FE field defines the FE approximation: 2_3_P1 = 2D, P1 on triangles.
field_1 = {
    'name' : 'temperature',
    'dim' : (1,1),
    'domain' : 'Omega',
    'bases' : {'Omega' : '2_3_P1'}
}

# Unknown and test functions (FE sense).
variables = {
    't' : ('unknown field', 'temperature', 0),
    's' : ('test field', 'temperature', 't'),
}

# Dirichlet boundary conditions.
ebcs = {
    't1' : ('Gamma_1', {'t.0' : 2.0}),
    't2' : ('Gamma_2', {'t.0' : -2.0}),
}

# Material coefficient c and source term value f.
material_1 = {
    'name' : 'coef',
    'mode' : 'here',
    'region' : 'Omega',
    'val' : 1.0,
}
material_2 = {
    'name' : 'source',
    'mode' : 'here',
    'region' : 'Omega_1',
    'val' : 10.0,
}

# Numerical quadrature and the equation.
integral_1 = {
    'name' : 'i1',
    'kind' : 'v',
    'quadrature' : 'gauss_o2_d2',
}

equations = {
    'Poisson' : """dw_laplace.i1.Omega( coef.val, s, t )
                 = dw_volume_lvf.i1.Omega_1( source.val, s )"""
}

# Solvers.
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
    'matrix'    : 'internal', # 'external' or 'internal'
    'problem'   : 'nonlinear', # 'nonlinear' or 'linear' (ignore i_max)
}

# FE assembling parameters: 'chunk_size' determines maximum number of elements
# to assemble in one C function call. Higher values mean faster assembling, but
# also more memory usage.
fe = {
    'chunk_size' : 10000
}

# Functions.
def select_circ( x, y, z, diameter ):
    """Select circular subdomain of a given diameter."""
    r = nm.sqrt( x**2 + y**2 )

    out = nm.where( r < diameter, 1, 0 )

    n = nm.where( out == 1 )[0].shape[0]
    if n <= 3:
        raise ValueError( 'too few nodes selected! (%d)' % n )

    return out

def vary_omega1_size( problem ):
    """Vary size of \Omega1. Saves also the regions into options['output_dir'].

    Input:
      problem: ProblemDefinition instance
    Return:
      a generator object:
      1. creates new (modified) problem
      2. yields the new (modified) problem and output container
      3. use the output container for some logging
      4. yields None (to signal next iteration to Application)
    """
    from sfepy.fem.problemDef import ProblemDefinition
    from sfepy.solvers.ts import get_print_info
    
    set_output_prefix( 'vary_omega1_size:' )

    diameters = nm.linspace( 0.1, 0.6, 7 ) + 0.001
    ofn_trunk, output_dir = problem.ofn_trunk, problem.output_dir
    join = os.path.join

    conf = problem.conf
    cr = conf.get_raw( 'regions' )
    n_digit, aux, suffix = get_print_info( len( diameters ) + 1 )
    d_format = suffix[1:-4]
    for ii, diameter in enumerate( diameters ):
        output( 'iteration %d: diameter %3.2f' % (ii, diameter) )

        cr['Omega_1'] = ('nodes by select_circ( x, y, z, %.5f )' % diameter, {})
        conf.edit( 'regions', cr )
        problem = ProblemDefinition.from_conf( conf )

        problem.save_regions( join( output_dir, ('regions_' + d_format) % ii ),
			      ['Omega_1'] )
        region = problem.domain.regions['Omega_1']
	if not region.has_cells_if_can():
	    raise ValueError( 'region %s has no cells!' % region.name )

        problem.ofn_trunk = ofn_trunk + '_' + (d_format % ii)

        out = []
        yield problem, out

        out_problem, vec_dp, data = out[-1]

        filename = join( output_dir,
                         ('log_%s.txt' % d_format) % ii )
	fd = open( filename, 'w' )
        log_item = '$r(\Omega_1)$: %f\n' % diameter
	fd.write( log_item )
	fd.write( 'solution:\n' )
	nm.savetxt( fd, vec_dp )
	fd.close()

        yield None

