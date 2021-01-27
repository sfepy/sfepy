from __future__ import print_function
from __future__ import absolute_import
import numpy as np
import six

try:
    import enthought.mayavi as mayavi

except ImportError:
    try:
        import mayavi

    except ImportError:
        mayavi = None

if mayavi:
    from mayavi import mlab
    from mayavi.core.source import Source
    from mayavi.core.filter import Filter
    import mayavi.core.dataset_manager as dataset_manager

from sfepy.base.base import basestr

def get_data_ranges(obj, return_only=False, use_names=None, filter_names=None):
    """Collect and print information on ranges of data in a dataset.

    Parameters
    ----------
    obj : a mayavi pipeline object
        The object to probe for data.

    return_only : boolean
        If True, do not print the information, just return it to the caller.

    use_names : list of strings
        Consider only data with names in the list.

    filter_names : list of strings
        Consider only data with names not in the list.

    Returns
    -------
    ranges : dict
        The requested data ranges.
    """
    if isinstance(obj, basestr):
        mlab.options.offscreen = True
        scene = mlab.figure(bgcolor=(1,1,1), fgcolor=(0, 0, 0), size=(1, 1))
        obj = mlab.pipeline.open(obj)

    if filter_names is None:
        filter_names = []

    source = obj
    while not (isinstance(source, Source)
               and not isinstance(source, Filter)):
        source = source.parent

    try:
        dm = dataset_manager.DatasetManager(dataset=source.outputs[0].output)

    except AttributeError:
        # Prior to Mayavi 4.6.2.
        try:
            dm = dataset_manager.DatasetManager(dataset=source.outputs[0])

        except:
            # Mayavi 4.7.1.
            dm = dataset_manager.DatasetManager(dataset=source.data)

    ranges = {}
    for attr in ['point_scalars', 'point_vectors', 'point_tensors',
                 'cell_scalars', 'cell_vectors', 'cell_tensors']:
        family, kind = attr.split('_')

        data = getattr(dm, attr)

        if use_names is None:
            names = list(data.keys())
        else:
            names = use_names

        if not return_only and (set(data.keys()).intersection(names)):
            print(family, kind)

        for key, arr in six.iteritems(data):
            if (key in filter_names) or (key not in names): continue

            shape = arr.shape
            if arr.ndim == 1:
                arr = arr[:,np.newaxis]
            norm = np.sqrt(np.sum(arr*arr, axis=1))

            ranges[key] = (family, kind, shape, np.min(arr), np.max(arr),
                           np.min(norm), np.max(norm))

            if not return_only:
                aux = (key,) + ranges[key][2:]
                print('"%s" %s range: %s %s l2 norm range: %s %s' % aux)

    return ranges
