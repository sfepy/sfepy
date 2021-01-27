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
from __future__ import absolute_import
from copy import copy

from sfepy.base.base import assert_, Struct
from sfepy.postprocess.utils import mlab
import six

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

def _get_scalars(color_name, color_kind, active):
    """
    Get scalars out of active data.
    """
    if color_kind == 'tensors':
        new_name = '|%s|' % color_name
        active = mlab.pipeline.set_active_attribute(active)
        active.point_tensors_name = color_name
        active = mlab.pipeline.extract_tensor_components(active)

    elif color_kind == 'vectors':
        new_name = '|%s|' % color_name
        active = mlab.pipeline.set_active_attribute(active)
        active.point_vectors_name = color_name
        active = mlab.pipeline.extract_vector_norm(active)

    elif color_kind == 'scalars':
        new_name = '%s' % color_name
        active = mlab.pipeline.set_active_attribute(active)
        active.point_scalars_name = color_name

    return new_name, active

def plot_displacements(source, ctp, bbox, position, family, kind, name,
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
        new_name, active = _get_scalars(color_name, color_kind, active)

    surf = mlab.pipeline.surface(active, opacity=opacity)
    surf.actor.actor.position = position

    return new_kind, new_name, active

def plot_warp_scalar(source, ctp, bbox, position, family, kind, name,
                     rel_scaling=1.0,
                     color_kind=None, color_name=None, opacity=1.0):
    """
    Show a 2D scalar field by displaying a colormap given by quantity
    `color_name` on the deformed mesh deformed by the scalar in the third
    dimension.

    Parameters
    ----------
    rel_scaling : float
        The relative scaling of scalar warp.
    color_kind : str, optional
        The kind of data determining the colormap.
    color_name : str, optional
        The name of data determining the colormap.
    opacity : float
        The surface plot opacity.
    """
    assert_(kind == 'scalars')

    if color_name is None:
        active = mlab.pipeline.set_active_attribute(source)

    else:
        active = mlab.pipeline.set_active_attribute(ctp)

    active.point_scalars_name = name
    active = mlab.pipeline.warp_scalar(active)
    active.filter.scale_factor = rel_scaling

    if color_name is None:
        new_kind = kind
        new_name = name

    else:
        new_kind = 'scalars'
        new_name, active = _get_scalars(color_name, color_kind, active)

    surf = mlab.pipeline.surface(active, opacity=opacity)
    surf.actor.actor.position = position

    return new_kind, new_name, active

def plot_velocity(source, ctp, bbox, position, family, kind, name,
                  seed='sphere', type='ribbon', integration_direction='both',
                  seed_scale=1.0, seed_resolution=20,
                  widget_enabled=True,
                  color_kind=None, color_name=None, opacity=1.0,
                  **kwargs):
    """
    Show velocity field by displaying streamlines and optionally a
    surface plot given by quantity `color_name`.

    Parameters
    ----------
    seed : one of ('sphere', 'point', 'line', 'plane')
        The streamline seed name.
    type : one of ('line', 'ribbon', 'tube')
        The streamline seed line type.
    integration_direction : one of ('forward', 'backward', 'both')
        The stream tracer integration direction.
    seed_scale : float
        The seed size scale.
    seed_resolution : int
        The number of seed points in a direction (depends on `seed`).
    widget_enabled : bool
        It True, the seed widget is visible and can be interacted with.
    color_kind : str, optional
        The kind of data determining the colormap.
    color_name : str, optional
        The name of data determining the colormap.
    opacity : float
        The surface plot opacity.
    **kwargs : dict
        Additional keyword arguments for attributes of
        `streamline.seed.widget`.
    """
    assert_(kind == 'vectors')

    active_v = mlab.pipeline.set_active_attribute(source)
    active_v.point_vectors_name = name

    active_n = mlab.pipeline.extract_vector_norm(active_v)

    s = mlab.pipeline.streamline
    streamline = s(active_n, seedtype=seed,
                   linetype=type,
                   seed_visible=True,
                   seed_scale=seed_scale,
                   integration_direction=integration_direction,
                   seed_resolution=seed_resolution)
    streamline.update_streamlines = True
    streamline.seed.widget.enabled = widget_enabled
    streamline.actor.actor.position = position

    for key, val in six.iteritems(kwargs):
        setattr(streamline.seed.widget, key, val)

    if color_name is None:
        active = active_n

    else:
        new_name, active = _get_scalars(color_name, color_kind, source)

    surf = mlab.pipeline.surface(active, opacity=opacity)
    surf.actor.actor.position = position

    if color_name is not None:
        mm = active.children[0]
        lm = mm.scalar_lut_manager
        scalar_bars = [['point', new_name, lm]]

        active_v.point_vectors_name = name # This is needed to have
                                           # colors by velocity!
        return kind, name, active_n, scalar_bars

    else:
        return kind, name, active_n

