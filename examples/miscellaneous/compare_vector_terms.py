r"""
Example without a physical relevance comparing new and old-style terms
with vector variables.

Find :math:`\ul{u}` (new style terms), :math:`\ul{r}` (old_style terms)
such that:

.. math::
    \int_{\Omega} D_{ijkl}\ e_{ij}(\ul{v}) e_{kl}(\ul{u})
    \int_{\Omega}\ul{v} \cdot \ul{u}
    = 0
    \;, \quad \forall \ul{v} \;,

    \int_{\Omega} D_{ijkl}\ e_{ij}(\ul{s}) e_{kl}(\ul{r})
    \int_{\Omega}\ul{s} \cdot \ul{r}
    = 0
    \;, \quad \forall \ul{s} \;.

The same values of :math:`\ul{u}`, :math:`\ul{r}` should be obtained.
"""
import os

from sfepy import data_dir
from sfepy.discrete.fem import MeshIO
from sfepy.mechanics.matcoefs import stiffness_from_lame

filename_mesh = data_dir + '/meshes/3d/cylinder.mesh'
#filename_mesh = data_dir + '/meshes/3d/cube_medium_hexa.mesh'

conf_dir = os.path.dirname(__file__)
io = MeshIO.any_from_filename(filename_mesh, prefix_dir=conf_dir)
bbox = io.read_bounding_box()

dd = bbox[1] - bbox[0]

xmin, ymin, zmin = bbox[0, :] + 1e-4 * dd
xmax, ymax, zmax = bbox[1, :] - 1e-4 * dd

def post_process(out, pb, state, extend=False):

    for vn in ['u', 'r']:
        try:
            val = pb.evaluate('dw_new_mass.2.Omega(%s, %s)'
                              % (vn, vn), verbose=False)
            print 'dw_new_mass', vn, val

            val = pb.evaluate('dw_new_lin_elastic.2.Omega(m.D, %s, %s)'
                              % (vn, vn), verbose=False)
            print 'dw_new_lin_elastic', vn, val

            val = pb.evaluate('dw_lin_elastic.2.Omega(m.D, %s, %s)'
                              % (vn, vn), verbose=False)
            print 'dw_lin_elastic', vn, val

        except:
            pass

    return out

options = {
    'nls' : 'newton',
    'ls' : 'ls',
    'post_process_hook' : 'post_process',
}

materials = {
    'm' : ({'D' : stiffness_from_lame(3, lam=0.0007, mu=0.0003),
            'one' : 1.0},),
}

regions = {
    'Omega' : 'all',
    'Gamma_Left' : ('vertices in (x < %f)' % xmin, 'facet'),
    'Gamma_Right' : ('vertices in (x > %f)' % xmax, 'facet'),
}

fields = {
    'displacements' : ('real', 'vector', 'Omega', 1),
}

variables = {
    'u' : ('unknown field', 'displacements', 0),
    'v' : ('test field',    'displacements', 'u'),
    'r' : ('unknown field', 'displacements', 1),
    's' : ('test field',    'displacements', 'r'),
}

ebcs = {
    'u1' : ('Gamma_Left', {'u.all' : 0.0}),
    'u2' : ('Gamma_Right', {'u.0' : 0.1 * (xmax - xmin),
                            'u.1' : 0.1 * (ymax - ymin)}),
    'r1' : ('Gamma_Left', {'r.all' : 0.0}),
    'r2' : ('Gamma_Right', {'r.0' : 0.1 * (xmax - xmin),
                            'r.1' : 0.1 * (ymax - ymin)}),
}

equations = {
    'new equation' :
    """dw_new_lin_elastic.2.Omega(m.D, v, u)
     + dw_new_mass.2.Omega(v, u)
     = 0""",
    'equation' :
    """dw_lin_elastic.2.Omega(m.D, s, r)
     + dw_volume_dot.2.Omega(m.one, s, r)
     = 0""",
}

solvers = {
    'ls' : ('ls.scipy_direct', {}),
    'newton' : ('nls.newton',
                {'i_max'      : 1,
                 'eps_a'      : 1e-10,
    }),
}
