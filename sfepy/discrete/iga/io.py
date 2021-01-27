"""
IO for NURBS and Bezier extraction data.
"""
from __future__ import absolute_import
import numpy as nm
import six
from six.moves import range
from sfepy.base.ioutils import HDF5ContextManager, enc, dec

def write_iga_data(filename, group, knots, degrees, control_points, weights,
                   cs, conn, bezier_control_points, bezier_weights, bezier_conn,
                   regions, name=None):
    """
    Write IGA-related data into a HDF5 file using pytables.

    filename: str or tables.File
        File to read the hdf5 mesh to.
    group: tables.group.Group, optional
        HDF5 file group to read the data from.
        If None, the root of file is used.

    Returns
    -------
    tuple
        Data for restoring IGA domain.
    """

    with HDF5ContextManager(filename, mode = 'w',
                            title='SfePy IGA data file') as fd:
        if group is None:
            group = fd.root

        if isinstance(degrees, int): degrees = [degrees]
        degrees = nm.asarray(degrees)

        nurbs = fd.create_group(group, 'nurbs', 'nurbs')

        fd.create_array(nurbs, 'dim', (int)(control_points.shape[1]), 'dim')
        fd.create_array(nurbs, 'tdim', len(degrees), 'tdim')
        for ii, kv in enumerate(knots):
            key = 'knots_%d' % ii
            fd.create_array(nurbs, key, kv, key)
        fd.create_array(nurbs, 'degrees', degrees, 'degrees')
        fd.create_array(nurbs, 'control_points', control_points,
                        'control_points')
        fd.create_array(nurbs, 'weights', weights, 'weights')

        bezier = fd.create_group(group, 'bezier', 'bezier')

        fd.create_array(bezier, 'bezier_control_points', bezier_control_points,
                        'bezier_control_points')
        fd.create_array(bezier, 'bezier_weights', bezier_weights,
                        'bezier_weights')
        for ii, op in enumerate(cs):
            key = 'extraction_%d' % ii
            fd.create_array(bezier, key, op, key)
        fd.create_array(bezier, 'global_connectivity', conn,
                        'global_connectivity')
        fd.create_array(bezier, 'bezier_connectivity', bezier_conn,
                        'bezier_connectivity')

        regs = fd.create_group(group, 'regions', 'regions')
        for key, val in six.iteritems(regions):
            fd.create_array(regs, key, val, key)

        if name is not None:
            fd.create_array( group, 'name', nm.array( enc(name)) )

def read_iga_data(filename, group=None):
    """
    Read IGA-related data from a HDF5 file using pytables.

    filename: str or tables.File
        File to read the hdf5 mesh to.
    group: tables.group.Group or None
        HDF5 file group to read the mesh from.
        If it's None, the root of file is used.

    Returns
    -------
    tuple
       Data for restoring IGA domain.
    """

    with HDF5ContextManager(filename, 'r') as fd:
        if group is None:
            group = fd.root

        nurbs = group.nurbs

        tdim = nurbs.tdim.read()

        knots = []
        for ii in range(tdim):
            name = 'knots_%d' % ii
            knots.append(nurbs._f_get_child(name).read())
        knots = tuple(knots)

        degrees = nurbs.degrees.read()
        control_points = nurbs.control_points.read()
        weights = nurbs.weights.read()

        bezier = group.bezier

        cs = []
        for ii in range(tdim):
            name = 'extraction_%d' % ii
            cs.append(bezier._f_get_child(name).read())

        conn = bezier.global_connectivity.read()
        bezier_control_points = bezier.bezier_control_points.read()
        bezier_weights = bezier.bezier_weights.read()
        bezier_conn = bezier.bezier_connectivity.read()

        regions = {}
        for region in group.regions:
            regions[region.name] = region.read()

        out = (knots, degrees, control_points, weights, cs, conn,
               bezier_control_points, bezier_weights, bezier_conn, regions)

        if hasattr(group, 'name'):
           out = out + ( dec(group.name.read().item()), )

        return out
