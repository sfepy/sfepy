try:
    from extmods import *
except (ImportError, AttributeError):
    from sfepy.base.base import output
    msg = 'sfepy extension modules are not compiled!\ntype "make"'
    raise ImportError( msg )

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
