# FEM stuff to move/split to common.
from __future__ import absolute_import
from sfepy.discrete.common.domain import Domain
from sfepy.discrete.common.region import Region
from sfepy.discrete.common.fields import Field

from .functions import Functions, Function
from .conditions import Conditions
from .variables import Variables, Variable, FieldVariable, DGFieldVariable, create_adof_conns
from .materials import Materials, Material
from .equations import Equations, Equation
from .integrals import Integrals, Integral
from .state import State
from .problem import Problem
from .evaluate import assemble_by_blocks
