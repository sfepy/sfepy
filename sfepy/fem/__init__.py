try:
    from extmods import *
except:
    from sfepy.base.base import output
    output( 'warning: sfepy extension modules are not compiled!' )
    output( 'type "make"' )


from mesh import Mesh
from domain import Domain
from fields import Fields, Field
from variables import Variables, Variable
from materials import Materials, Material
from equations import Equations, Equation
from integrals import Integrals, Integral
from problemDef import ProblemDefinition
from sfepy.fem.meshio import MeshIO
from evaluate import eval_term_op
