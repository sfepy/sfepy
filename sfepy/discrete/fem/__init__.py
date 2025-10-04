try:
    from .extmods import *
except (ImportError, AttributeError):
    print('sfepy extension modules may not be compiled!\nTry typing "make".')
    raise

from .mesh import Mesh
from .domain import FEDomain
from .fields_base import Field
from sfepy.discrete.fem.meshio import MeshIO
from .utils import extend_cell_data
