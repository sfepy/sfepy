import numpy as nm
from linear_elastic import \
     filename_mesh, materials, regions, fields, ebcs, \
     integrals, solvers, fe

options = {
    'ts' : 'ts',
    'save_steps' : -1,
}

variables = {
    'u' : ('unknown field', 'displacement', 0, 'previous'),
    'v' : ('test field', 'displacement', 'u'),
}

# Put density to 'solid'.
materials['solid'][0].update({'rho' : 1000.0})

# Moving the PerturbedSurface region.
ebcs['PerturbedSurface'][1].update({'u.0' : 'ebc_sin'})

def ebc_sin(ts, coors, bc=None):
    val = 0.01 * nm.sin(2.0*nm.pi*ts.nt)
    return nm.tile(val, (coors.shape[0],))

functions = {
    'ebc_sin' : (ebc_sin,),
}

equations = {
    'balance_of_forces in time' :
    """dw_mass_vector.i1.Omega( solid.rho, v, du/dt )
     + dw_lin_elastic_iso.i1.Omega( solid.lam, solid.mu, v, u ) = 0""",
}

solvers.update({
    'ts' : ('ts.simple',
            {'t0' : 0.0,
             't1' : 1.0,
             'dt' : None,
             'n_step' : 101
             }),
})

# Pre-assemble and factorize the matrix prior to time-stepping.
newton = solvers['newton']
newton[1].update({'problem' : 'linear'})

ls = solvers['ls']
ls[1].update({'presolve' : True})
