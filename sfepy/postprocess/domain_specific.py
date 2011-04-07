"""
Domain-specific plot functions.

All the plot functions accept the following parameters:

- `source` : Mayavi source
- `ctp` : Mayavi cell-to-point filter
- `position` : `(x, y, z)`
- `family` : 'point' or 'cell'
- `kind` : 'scalars', 'vectors' or 'tensors'
- `name` : name of a variable

All the plot functions return:
- `kind` : 'scalars', 'vectors' or 'tensors'
- `name` : name of a variable
- `active` : Mayavi module
"""
from copy import copy

from sfepy.base.base import assert_, Struct
from sfepy.postprocess.utils import mlab

class DomainSpecificPlot(Struct):
    """
    Class holding domain-specific plot function and its parameters.
    """
    def __init__(self, fun_name, args):
        Struct.__init__(self, fun_name=fun_name)

        ii = fun_name.rfind('.')
        if ii >= 0:
            # Function in a user module.
            mod_name, fun_name = fun_name[:ii], fun_name[ii+1:]
            mod = __import__(mod_name, fromlist=[fun_name])
            self.fun = getattr(mod, fun_name)

        else:
            self.fun = globals()[self.fun_name]

        kwargs = {}
        for arg in args:
            key, val = arg.split('=')
            kwargs[key] = eval(val)

        self.kwargs = kwargs

    def __call__(self, *args, **kwargs):
        _kwargs = copy(kwargs)
        _kwargs.update(self.kwargs)

        return self.fun(*args, **_kwargs)

def plot_displacements(source, ctp, position, family, kind, name,
                       rel_scaling=1.0,
                       color_kind=None, color_name=None, opacity=1.0):
    """
    Show displacements by displaying a colormap given by quantity
    `color_name` on the deformed mesh.

    Parameters
    ----------
    rel_scaling : float
        The relative scaling of displacements.
    color_kind : str, optional
        The kind of data determining the colormap.
    color_name : str, optional
        The name of data determining the colormap.
    opacity : float
        The surface plot opacity.
    """
    assert_(kind == 'vectors')

    if color_name is None:
        active = mlab.pipeline.set_active_attribute(source)

    else:
        active = mlab.pipeline.set_active_attribute(ctp)

    active.point_vectors_name = name
    active = mlab.pipeline.warp_vector(active)
    active.filter.scale_factor = rel_scaling

    if color_name is None:
        new_kind = kind
        new_name = name

        active = mlab.pipeline.extract_vector_norm(active)

    else:
        new_kind = 'scalars'

        if color_kind == 'tensors':
            new_name = '|%s|' % color_name
            active = mlab.pipeline.set_active_attribute(active)
            active.point_tensors_name = color_name
            active = mlab.pipeline.extract_tensor_components(active)

        elif color_kind == 'vectors':
            new_name = '|%s|' % color_name
            active = mlab.pipeline.set_active_attribute(active)
            active.point_vectors_name = color_name
            active = mlab.pipeline.extract_tensor_components(active)

        elif color_kind == 'scalars':
            new_name = '%s' % color_name
            active = mlab.pipeline.set_active_attribute(active)
            active.point_scalars_name = color_name

    surf = mlab.pipeline.surface(active, opacity=opacity)
    surf.actor.actor.position = position

    return new_kind, new_name, active
