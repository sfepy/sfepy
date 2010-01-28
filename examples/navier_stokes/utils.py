##
# Functions.
import numpy as nm

def coors_in_cylinder( coors, centre, axis, radius, length ):
    """
    Select coordinates in a cylindre.
    """
    vec = coors.T - centre
##     ort_axis = nm.cross( axis, axis + axis[[2,1,0]] * [2,3,4] )
##     norm = nm.sqrt( nm.dot( ort_axis, ort_axis ) )
##     ort_axis = ort_axis / norm
    
    drv = nm.cross( axis, vec, axisb = 0 )
    dr = nm.sqrt( nm.sum( drv * drv, 1 ) )
    dl = nm.dot( axis, vec )
    out = nm.where( (dl >= 0.0) & (dl <= length) & (dr <= radius))[0]

    return out

# last revision: 01.08.2007
def cinc_cylinder( coors, mode ):
    axis = nm.array( [1, 0, 0], nm.float64 )
    if mode == 0: # In
        centre = nm.array( [-0.00001, 0.0, 0.0], nm.float64 ).reshape( (3,1) )
        radius = 0.019
        length = 0.00002
    elif mode == 1: # Out
        centre = nm.array( [0.09999, 0.0, 0.0], nm.float64 ).reshape( (3,1) )
        radius = 0.019
        length = 0.00002
    else:
        centre = nm.array( [0.05, 0.0, 0.0], nm.float64 ).reshape( (3,1) )
        radius = 0.012
        length = 0.04

    return coors_in_cylinder( coors, centre, axis, radius, length )

def cinc_elbow2( coors, mode ):
    if mode == 0: # In
        centre = nm.array( [0.0, -0.00001, 0.0], nm.float64 ).reshape( (3,1) )
    else: # Out
        centre = nm.array( [0.2, -0.00001, 0.0], nm.float64 ).reshape( (3,1) )
    
    axis = nm.array( [0, 1, 0], nm.float64 )
    radius = 0.029
    length = 0.00002

    return coors_in_cylinder( coors, centre, axis, radius, length )
