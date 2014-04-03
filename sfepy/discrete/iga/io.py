"""
IO for NURBS and Bezier extraction data.
"""
import numpy as nm
import tables as pt

def write_iga_data(filename, knots, degrees, control_points, weights, cs, conn,
                   bezier_control_points, bezier_weights, bezier_conn, regions):
    """
    Write IGA-related data into a HDF5 file using pytables.
    """
    if isinstance(degrees, int): degrees = [degrees]
    degrees = nm.asarray(degrees)

    fd = pt.openFile(filename, mode='w', title='SfePy IGA data file')

    nurbs = fd.createGroup('/', 'nurbs', 'nurbs')

    fd.createArray(nurbs, 'dim', control_points.shape[1], 'dim')
    fd.createArray(nurbs, 'tdim', len(degrees), 'tdim')
    for ii, kv in enumerate(knots):
        name = 'knots_%d' % ii
        fd.createArray(nurbs, name, kv, name)
    fd.createArray(nurbs, 'degrees', degrees, 'degrees')
    fd.createArray(nurbs, 'control_points', control_points, 'control_points')
    fd.createArray(nurbs, 'weights', weights, 'weights')

    bezier = fd.createGroup('/', 'bezier', 'bezier')

    fd.createArray(bezier, 'bezier_control_points', bezier_control_points,
                   'bezier_control_points')
    fd.createArray(bezier, 'bezier_weights', bezier_weights, 'bezier_weights')
    for ii, op in enumerate(cs):
        name = 'extraction_%d' % ii
        fd.createArray(bezier, name, op, name)
    fd.createArray(bezier, 'global_connectivity', conn, 'global_connectivity')
    fd.createArray(bezier, 'bezier_connectivity', bezier_conn,
                   'bezier_connectivity')

    regs = fd.createGroup('/', 'regions', 'regions')
    for key, val in regions.iteritems():
        fd.createArray(regs, key, val, key)

    fd.close()

def read_iga_data(filename):
    """
    Read IGA-related data from a HDF5 file using pytables.
    """
    fd = pt.openFile(filename, mode='r')

    nurbs = fd.root.nurbs

    tdim = nurbs.tdim.read()

    knots = []
    for ii in xrange(tdim):
        name = 'knots_%d' % ii
        knots.append(nurbs._f_getChild(name).read())
    knots = tuple(knots)

    degrees = nurbs.degrees.read()
    control_points = nurbs.control_points.read()
    weights = nurbs.weights.read()

    bezier = fd.root.bezier

    cs = []
    for ii in xrange(tdim):
        name = 'extraction_%d' % ii
        cs.append(bezier._f_getChild(name).read())

    conn = bezier.global_connectivity.read()
    bezier_control_points = bezier.bezier_control_points.read()
    bezier_weights = bezier.bezier_weights.read()
    bezier_conn = bezier.bezier_connectivity.read()

    regions = {}
    for region in fd.root.regions:
        regions[region.name] = region.read()

    fd.close()

    return (knots, degrees, control_points, weights, cs, conn,
            bezier_control_points, bezier_weights, bezier_conn, regions)
