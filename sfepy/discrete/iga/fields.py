"""
Fields for isogeometric analysis.
"""
import numpy as nm

from sfepy.base.base import Struct
from sfepy.discrete.common.fields import parse_shape, Field

class IGField(Field):
    """
    Bezier extraction based NURBS field for isogeometric analysis.
    """
    family_name = 'volume_H1_iga'

    def __init__(self, name, dtype, shape, region, **kwargs):
        """
        Create a Bezier element isogeometric analysis field.

        Parameters
        ----------
        name : str
            The field name.
        dtype : numpy.dtype
            The field data type: float64 or complex128.
        shape : int/tuple/str
            The field shape: 1 or (1,) or 'scalar', space dimension (2, or (2,)
            or 3 or (3,)) or 'vector', or a tuple. The field shape determines
            the shape of the FE base functions and is related to the number of
            components of variables and to the DOF per node count, depending
            on the field kind.
        region : Region
            The region where the field is defined.
        **kwargs : dict
            Additional keyword arguments.
        """
        shape = parse_shape(shape, region.domain.shape.dim)
        Struct.__init__(self, name=name, dtype=dtype, shape=shape,
                        region=region)

        self.domain = self.region.domain
        self.nurbs = self.domain.nurbs

        self._setup_kind()

        self.n_components = nm.prod(self.shape)
        self.val_shape = self.shape
        self.n_nod = self.nurbs.weights.shape[0]

        self.mappings = {}

        self.igs = self.region.igs

    def get_true_order(self):
        return nm.prod(self.nurbs.degrees)

    def is_higher_order(self):
        """
        Return True, if the field's approximation order is greater than one.
        """
        return (self.nurbs.degrees > 1).any()
