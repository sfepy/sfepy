##
# Functions.
from __future__ import absolute_import
import numpy as nm

from sfepy.linalg import get_coors_in_tube

# last revision: 01.08.2007
def cinc_cylinder(coors, mode):
    axis = nm.array([1, 0, 0], nm.float64)
    if mode == 0: # In
        centre = nm.array([-0.00001, 0.0, 0.0], nm.float64)
        radius = 0.019
        length = 0.00002
    elif mode == 1: # Out
        centre = nm.array([0.09999, 0.0, 0.0], nm.float64)
        radius = 0.019
        length = 0.00002
    else:
        centre = nm.array([0.05, 0.0, 0.0], nm.float64)
        radius = 0.012
        length = 0.04

    return get_coors_in_tube(coors, centre, axis, -1.0, radius, length)

def cinc_elbow2(coors, mode):
    if mode == 0: # In
        centre = nm.array([0.0, -0.00001, 0.0], nm.float64)
    else: # Out
        centre = nm.array([0.2, -0.00001, 0.0], nm.float64)
    
    axis = nm.array([0, 1, 0], nm.float64)
    radius = 0.029
    length = 0.00002

    return get_coors_in_tube(coors, centre, axis, -1.0, radius, length)
