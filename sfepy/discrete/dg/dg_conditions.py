from __future__ import absolute_import
import numpy as nm

from sfepy.base.base import basestr, Container, Struct
from sfepy.discrete.functions import Function

from sfepy.discrete.conditions import Condition, PeriodicBC, EssentialBC
import six


class DGPeriodicBC(PeriodicBC):
    ...

class DGEssentialBC(EssentialBC):
    ...