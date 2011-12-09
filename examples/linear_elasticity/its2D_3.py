r"""
Diametrically point loaded 2-D disk with nodal stress calculation. See
:ref:`sec-primer`.

Find :math:`\ul{u}` such that:

.. math::
    \int_{\Omega} D_{ijkl}\ e_{ij}(\ul{v}) e_{kl}(\ul{u})
    = 0
    \;, \quad \forall \ul{v} \;,

where

.. math::
    D_{ijkl} = \mu (\delta_{ik} \delta_{jl}+\delta_{il} \delta_{jk}) +
    \lambda \ \delta_{ij} \delta_{kl}
    \;.
"""
from its2D_1 import *

from sfepy.mechanics.matcoefs import stiffness_tensor_youngpoisson
from sfepy.fem.geometry_element import geometry_data
from sfepy.fem import Field,FieldVariable
import numpy as nm

gdata = geometry_data['2_3']
nc = len(gdata.coors)

def nodal_stress(out, pb, state, extend=False):
    '''
    Calculate stresses at nodal points.
    '''

    # Calc point load
    u = state.vec
    pb.remove_bcs()
    f = pb.evaluator.eval_residual(u)
    f.shape = (pb.domain.mesh.n_nod,2)
    P = 2. * f[2][1]

    # Calc nodal stress
    pb.time_update()
    ev = pb.evaluate

    stress = ev('dq_cauchy_stress.ivn.Omega(Asphalt.D,u)', mode='el_avg')
    sfield = Field('stress', nm.float64, (3,), pb.domain.regions["Omega"])
    svar = FieldVariable('sigma', 'parameter', sfield, 3,
                         primary_var_name='(set-to-None)')
    svar.data_from_qp(stress, pb.integrals['ivn'])

    print '\n=================================================================='
    print 'Load to give 1 mm displacement = %s Newton ' % round(-P,3)
    print '\nAnalytical solution approximation'
    print '================================='
    print 'Horizontal tensile stress = %s MPa/mm' % round(-2.*P/(nm.pi*150.),3) 
    print 'Vertical compressive stress = %s MPa/mm' % round(-6.*P/(nm.pi*150.),3)
    print '\nFEM solution'
    print '============'
    print 'Horizontal tensile stress = %s MPa/mm' % round(svar()[0][0],3)
    print 'Vertical compressive stress = %s MPa/mm' % -round(svar()[0][1],3)
    print '=================================================================='
    return out

asphalt = materials['Asphalt'][0]
asphalt.update({'D' : stiffness_tensor_youngpoisson(2, young, poisson)})
options.update({'post_process_hook' : 'nodal_stress',})

integrals = {
    'ivn' : ('v', 'custom', gdata.coors, [gdata.volume / nc] * nc),
}
