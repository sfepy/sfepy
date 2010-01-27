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
def cinc_simple( coors, mode ):
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

cinc_simple2 = cinc_simple
cinc_simple3 = cinc_simple

def cinc_klikatak( coors ):
    raise NotImplementedError

    return coors_in_cylinder( coors, centre, axis, radius, length )

def cinc_pul_klikatak2( coors, mode ):
    if mode == 0: # In
        centre = nm.array( [0.0, -0.00001, 0.0], nm.float64 ).reshape( (3,1) )
    else: # Out
        centre = nm.array( [0.2, -0.00001, 0.0], nm.float64 ).reshape( (3,1) )
    
    axis = nm.array( [0, 1, 0], nm.float64 )
    radius = 0.029
    length = 0.00002

    return coors_in_cylinder( coors, centre, axis, radius, length )

##
# 21.03.2007, c
# last revision: 21.03.2007
def cinc_new_kroucenak_c1( coors, mode ):
    if mode == 0: # In
        centre = nm.array( [-35.789001e-3, -37.777100e-3, -0.00001e-3],
                           nm.float64 ).reshape( (3,1) )
    else: # Out
        centre = nm.array( [13.083601e-3, -12.341600e-3, -0.00001e-3],
                           nm.float64 ).reshape( (3,1) )
    
    axis = nm.array( [0, 0, 1], nm.float64 )
    radius = 4.5e-3
    length = 0.00002e-3

    return coors_in_cylinder( coors, centre, axis, radius, length )

##
# 22.03.2007, c
# last revision: 22.03.2007
def cinc_new_kroucenak_c2( coors, mode ):
    if mode == 0: # In
        centre = nm.array( [-0.036469, -0.013616, -0.00001e-3],
                           nm.float64 ).reshape( (3,1) )
    else: # Out
        centre = nm.array( [0.013664, -0.036469, -0.00001e-3],
                           nm.float64 ).reshape( (3,1) )
    
    axis = nm.array( [0, 0, 1], nm.float64 )
    radius = 4.5e-3
    length = 0.00002e-3

    return coors_in_cylinder( coors, centre, axis, radius, length )

##
# 26.03.2007, c
# last revision: 26.03.2007
def cinc_new_kroucenak_c3( coors, mode ):
    if mode == 0: # In
        centre = nm.array( [-0.038927, -0.025294, -0.00001e-3],
                           nm.float64 ).reshape( (3,1) )
    else: # Out
        centre = nm.array( [0.016166, -0.024818, -0.00001e-3],
                           nm.float64 ).reshape( (3,1) )
    
    axis = nm.array( [0, 0, 1], nm.float64 )
    radius = 4.5e-3
    length = 0.00002e-3

    return coors_in_cylinder( coors, centre, axis, radius, length )
