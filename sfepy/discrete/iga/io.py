"""
IO for NURBS and Bezier extraction data.
"""
from __future__ import absolute_import
import numpy as nm
import tables as pt
import six
from six.moves import range

def write_iga_data(filename, knots, degrees, control_points, weights, cs, conn,
                   bezier_control_points, bezier_weights, bezier_conn, regions):
    """
    Write IGA-related data into a HDF5 file using pytables.
    """
    if isinstance(degrees, int): degrees = [degrees]
    degrees = nm.asarray(degrees)

    fd = pt.open_file(filename, mode='w', title='SfePy IGA data file')

    nurbs = fd.create_group('/', 'nurbs', 'nurbs')

    fd.create_array(nurbs, 'dim', control_points.shape[1], 'dim')
    fd.create_array(nurbs, 'tdim', len(degrees), 'tdim')
    for ii, kv in enumerate(knots):
        name = 'knots_%d' % ii
        fd.create_array(nurbs, name, kv, name)
    fd.create_array(nurbs, 'degrees', degrees, 'degrees')
    fd.create_array(nurbs, 'control_points', control_points, 'control_points')
    fd.create_array(nurbs, 'weights', weights, 'weights')

    bezier = fd.create_group('/', 'bezier', 'bezier')

    fd.create_array(bezier, 'bezier_control_points', bezier_control_points,
                   'bezier_control_points')
    fd.create_array(bezier, 'bezier_weights', bezier_weights, 'bezier_weights')
    for ii, op in enumerate(cs):
        name = 'extraction_%d' % ii
        fd.create_array(bezier, name, op, name)
    fd.create_array(bezier, 'global_connectivity', conn, 'global_connectivity')
    fd.create_array(bezier, 'bezier_connectivity', bezier_conn,
                   'bezier_connectivity')

    regs = fd.create_group('/', 'regions', 'regions')
    for key, val in six.iteritems(regions):
        fd.create_array(regs, key, val, key)

    fd.close()

def read_iga_data(filename):
    """
    Read IGA-related data from a HDF5 file using pytables.
    """
    fd = pt.open_file(filename, mode='r')

    nurbs = fd.root.nurbs

    tdim = nurbs.tdim.read()

    knots = []
    for ii in range(tdim):
        name = 'knots_%d' % ii
        knots.append(nurbs._f_get_child(name).read())
    knots = tuple(knots)

    degrees = nurbs.degrees.read()
    control_points = nurbs.control_points.read()
    weights = nurbs.weights.read()

    bezier = fd.root.bezier

    cs = []
    for ii in range(tdim):
        name = 'extraction_%d' % ii
        cs.append(bezier._f_get_child(name).read())

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
