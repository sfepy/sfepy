import numpy as np

try:
    from enthought.mayavi import mlab
except ImportError:
    mlab = None

if mlab:
    from enthought.mayavi.core.source import Source
    from enthought.mayavi.core.filter import Filter
    import dataset_manager

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
    if isinstance(obj, str):
        mlab.options.offscreen = True
        scene = mlab.figure(bgcolor=(1,1,1), fgcolor=(0, 0, 0), size=(1, 1))
        obj = mlab.pipeline.open(obj)
    
    if filter_names is None:
        filter_names = []

    source = obj
    while not (isinstance(source, Source)
               and not isinstance(source, Filter)):
        source = source.parent

    dm = dataset_manager.DatasetManager(dataset=source.outputs[0])

    ranges = {}
    for attr in ['point_scalars', 'point_vectors', 'point_tensors',
                 'cell_scalars', 'cell_vectors', 'cell_tensors']:
        family, kind = attr.split('_')

        data = getattr(dm, attr)

        if use_names is None:
            names = data.keys()
        else:
            names = use_names

        if not return_only and (set(data.keys()).intersection(names)):
            print family, kind

        for key, arr in data.iteritems():
            if (key in filter_names) or (key not in names): continue

            shape = arr.shape
            if arr.ndim == 1:
                arr = arr[:,np.newaxis]
            norm = np.sqrt(np.sum(arr*arr, axis=1))

            ranges[key] = (family, kind, shape, np.min(arr), np.max(arr),
                           np.min(norm), np.max(norm))

            if not return_only:
                aux = (key,) + ranges[key][2:]
                print '  "%s" %s range: %s %s l2_norm_range: %s %s' % aux

    return ranges
