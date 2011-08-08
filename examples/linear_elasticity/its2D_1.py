# Diametrically point loaded 2-D disk.

from sfepy.mechanics.matcoefs import youngpoisson_to_lame
from sfepy import data_dir

filename_mesh = data_dir + '/meshes/2d/its2D.mesh'

output_dir = '.' # set this to a valid directory you have write access to

young = 2000.0 # Young's modulus [MPa]
poisson = 0.4  # Poisson's ratio

options = {
    'output_dir' : output_dir,
}

regions = {
    'Omega' : ('all', {}),
    'Left' : ('nodes in (x < 0.001)', {}),
    'Bottom' : ('nodes in (y < 0.001)', {}),
    'Top' : ('node 2', {}),
}

materials = {
    'Asphalt' : ({
        'lam' : youngpoisson_to_lame(young, poisson)[0],
        'mu' : youngpoisson_to_lame(young, poisson)[1],
    },),
}

fields = {
    'displacement': ('real', 'vector', 'Omega', 1),
}

equations = {
   'balance_of_forces' :
   """dw_lin_elastic_iso.2.Omega(Asphalt.lam, Asphalt.mu, v, u ) = 0""",
}

variables = {
    'u' : ('unknown field', 'displacement', 0),
    'v' : ('test field', 'displacement', 'u'),
}

ebcs = {
    'XSym' : ('Bottom', {'u.1' : 0.0}),
    'YSym' : ('Left', {'u.0' : 0.0}),
    'Load' : ('Top', {'u.0' : 0.0, 'u.1' : -1.0}),
}

solvers = {
    'ls' : ('ls.scipy_direct', {}),
    'newton' : ('nls.newton', {
        'i_max' : 1,
        'eps_a' : 1e-6,
        'problem' : 'nonlinear'
    }),
}

