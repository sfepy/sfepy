import os
from sfepy import top_dir

# Homogenization - get material properties
def get_homog( ts, coor, mode, region, ig,
               filename = None, coefs_filename=None,
               regenerate = False ):
     
     import tables as pt
     import numpy as nm
     
     if not regenerate:
         if os.path.exists( coefs_filename ):
             if not pt.isHDF5File( coefs_filename ):
                 regenerate = True
         else:
             regenerate = True
   
     if regenerate:
         from sfepy.base.base import Struct, output
         from sfepy.base.conf import ProblemConf, get_standard_keywords
         from homogen import HomogenizationApp
         
         output.prefix = 'micro:'

         required, other = get_standard_keywords()
         required.remove( 'equations' )
         
         conf = ProblemConf.from_file( filename, required, other )
         options = Struct( output_filename_trunk = None )
   
         app = HomogenizationApp( conf, options, 'micro:' )
         coefs = app()
         
         coefs.to_file_hdf5( coefs_filename )
         
         output.prefix = 'macro:'
     else:
         from sfepy.homogenization.coefficients import Coefficients
         coefs = Coefficients.from_file_hdf5( coefs_filename )

     out = {}
     if mode == 'qp':
         out['D'] = nm.tile( coefs.D, (coor.shape[0], 1, 1) )
         
     return  out


functions = {
    'get_homog' : (lambda ts, coors, mode = None, region = None, ig = None:
                   get_homog( ts, coors, mode, region, ig,
                              filename = 'examples/homogenization/linear_homogenization.py',
                              coefs_filename = './output/coefs_le.h5',
                              regenerate = False ),),
}

filename_mesh = top_dir + '/meshes/3d/cylinder.mesh'

regions = {
    'Omega' : ('all', {}),
    'Left' : ('nodes in (x < 0.001)', {}),
    'Right' : ('nodes in (x > 0.099)', {}),
}

materials = {
    'solid' : ('Omega', None, 'get_homog'),
}

fields = {
    '3_displacement': ((3,1), 'real', 'Omega', {'Omega' : '3_4_P1'}),
}

integrals = {
    'i1' : ('v', 'gauss_o1_d3'),
}

variables = {
    'u' : ('unknown field', '3_displacement', 0),
    'v' : ('test field', '3_displacement', 'u'),
}

ebcs = {
    'Fixed' : ('Left', {'u.all' : 0.0}),
    'PerturbedSurface' : ('Right', {'u.0' : 0.02, 'u.1' : 0.0, 'u.2' : 0.0}),
}

equations = {
    'balance_of_forces' :
    """dw_lin_elastic.i1.Omega( solid.D, v, u ) = 0""",
}

solvers = {
    'ls' : ('ls.scipy_direct', {}),
    'newton' : ('nls.newton',
                { 'i_max'      : 1,
                  'eps_a'      : 1e-8,
                  'problem'   : 'nonlinear'}),
}

fe = {
    'chunk_size' : 10000
}

options = {
    'nls' : 'newton',
    'ls' : 'ls',
    'output_dir' : 'output',
}
