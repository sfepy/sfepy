from __future__ import absolute_import
import numpy as nm

from sfepy.base.base import basestr, Container, Struct
from sfepy.discrete.functions import Function

from sfepy.discrete.conditions import Condition, PeriodicBC, EssentialBC
import six


class DGPeriodicBC(PeriodicBC):
    """
    This class is empty, it serves the same purpose
    as PeriodicBC, and is created only for branching in
    dof_info.py
    """
    pass

class DGEssentialBC(EssentialBC):
    """
    This class is empty, it serves the same purpose
    as EssentialBC, and is created only for branching in
    dof_info.py
    """

    def __init__(self, *args, diff=0, **kwargs):
        EssentialBC.__init__(self, *args, **kwargs)
        self.diff = diff
