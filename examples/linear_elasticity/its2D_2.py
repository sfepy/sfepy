from its2D_1 import *

from sfepy.mechanics.matcoefs import stiffness_tensor_youngpoisson

def stress_strain(out, pb, state, extend=False):
    """
    Calculate and output strain and stress for given displacements.
    """
    from sfepy.base.base import Struct

    ev = pb.evaluate
    strain = ev('de_cauchy_strain.2.Omega(u)', mode='el_avg')
    stress = ev('de_cauchy_stress.2.Omega(Asphalt.D, u)', mode='el_avg')

    out['cauchy_strain'] = Struct(name='output_data', mode='cell',
                                  data=strain, dofs=None)
    out['cauchy_stress'] = Struct(name='output_data', mode='cell',
                                  data=stress, dofs=None)

    return out

asphalt = materials['Asphalt'][0]
asphalt.update({'D' : stiffness_tensor_youngpoisson(2, young, poisson)})
options.update({'post_process_hook' : 'stress_strain',})
